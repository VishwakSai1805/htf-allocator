/* ============================================================================
 *  Hybrid Threshold Fit (HTF) — a dynamic memory allocator in C
 *  ---------------------------------------------------------------------------
 *  OS course project.  Single file, terminal only, no external libraries.
 *
 *  WHAT THIS PROGRAM IS
 *  --------------------
 *  It sets aside one big static array (the "heap") and implements malloc/free
 *  on top of it, six different ways:
 *
 *      1. First Fit   - take the first free hole that is big enough
 *      2. Best  Fit   - take the smallest hole that is big enough
 *      3. Worst Fit   - take the largest hole
 *      4. Next  Fit   - like First Fit, but resume searching where we stopped
 *      5. Buddy       - round the request up to a power of two, split/merge
 *      6. HTF         - MY ALGORITHM: small requests come from fixed size-class
 *                       bins in O(1); large requests use Best Fit with
 *                       coalescing.  The split point is the threshold T.
 *
 *  It then runs the SAME workload through all six and prints how much memory
 *  each one wasted, how much searching it did, and how many requests failed.
 *
 *  WHY IT MATTERS
 *  --------------
 *  Repeated allocate/free leaves the heap full of scattered holes.  Eventually
 *  a large request fails even though plenty of memory is free - it is just in
 *  pieces too small to use.  That is fragmentation.  Real systems split small
 *  and large allocations to control it: Linux uses its buddy page allocator for
 *  big chunks and slab/slub for small objects.  This program asks where that
 *  split should be placed, and proves the answer with measurements.
 *
 *  THE BLOCK LAYOUT (used by fits 1-4 and by HTF)
 *  ----------------------------------------------
 *  Every block carries a header AND a footer ("boundary tags", Knuth).  The
 *  footer is what makes backward coalescing possible: when freeing a block we
 *  can look at the bytes immediately before our header and learn how big the
 *  previous block was, so we can merge with it.
 *
 *      offset ->  [ Header | payload .......... | Footer ]
 *                   24 B      size bytes           8 B
 *
 *  Walking forward: next_offset = offset + HDR + size + FTR.
 *  Walking back:    read the footer just before our header to get prev size.
 *
 *  DELIBERATE SCOPE DECISIONS (state these in the viva)
 *  ----------------------------------------------------
 *  - The heap has a FIXED size.  This is an experimental control, not a
 *    shortcut: an expandable heap would hide fragmentation, which is the exact
 *    effect being measured.
 *  - No thread safety.  Locking is orthogonal to the allocation policy; adding
 *    it would change none of the results.
 *
 *  Build:  gcc -O2 -Wall -Wextra -o htf htf_allocator.c
 *  Run:    ./htf              (prints every table)
 *          ./htf --csv        (also writes results.csv)
 * ==========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <assert.h>

/* ---------------------------------------------------------------------------
 * Configuration
 * -------------------------------------------------------------------------*/

#define HEAP_SIZE      (1u << 20)   /* 1 MB simulated heap                   */
#define ALIGNMENT      8u           /* every payload is 8-byte aligned       */
#define MIN_PAYLOAD    16u          /* no payload smaller than this          */

/* Buddy allocator orders: smallest block 32 B (2^5), largest 1 MB (2^20).   */
#define BUDDY_MIN_ORDER 5
#define BUDDY_MAX_ORDER 20

/* Size classes for the HTF bins.  Deliberately FINE-GRAINED, unlike buddy's
 * powers of two - that is precisely how HTF keeps internal waste small for
 * small objects.  A 33-byte request takes 48 bytes here, but 64 under buddy. */
#define NUM_CLASSES 10
static const size_t class_size[NUM_CLASSES] = {
    8, 16, 24, 32, 48, 64, 96, 128, 192, 256
};

/* ---------------------------------------------------------------------------
 * Policies
 * -------------------------------------------------------------------------*/

typedef enum {
    FIRST_FIT = 0, BEST_FIT, WORST_FIT, NEXT_FIT, BUDDY, HTF
} Policy;

static const char *policy_name[] = {
    "First Fit", "Best Fit", "Worst Fit", "Next Fit", "Buddy", "HTF"
};

/* ---------------------------------------------------------------------------
 * The heap itself, plus block tags
 * -------------------------------------------------------------------------*/

static unsigned char heap[HEAP_SIZE];

typedef struct {
    size_t size;        /* payload bytes, excluding the tags                */
    size_t requested;   /* what the caller actually asked for  -> internal
                         * fragmentation = size - requested                 */
    int    is_free;
    int    bin_class;   /* -1 = ordinary block; >=0 = owned by that HTF bin */
} Header;

typedef struct { size_t size; } Footer;

#define HDR ((size_t)sizeof(Header))
#define FTR ((size_t)sizeof(Footer))

static Header *hdr_at(size_t off)      { return (Header *)(heap + off); }
static void   *payload_at(size_t off)  { return heap + off + HDR; }
static size_t  off_of_payload(void *p) { return (size_t)((unsigned char *)p - heap) - HDR; }

static Footer *ftr_at(size_t off) {
    Header *h = hdr_at(off);
    return (Footer *)(heap + off + HDR + h->size);
}
static size_t block_bytes(size_t off)  { return HDR + hdr_at(off)->size + FTR; }
static size_t next_off(size_t off)     { return off + block_bytes(off); }

static size_t align_up(size_t n) { return (n + (ALIGNMENT - 1)) & ~(size_t)(ALIGNMENT - 1); }

/* ---------------------------------------------------------------------------
 * Global allocator state and counters
 * -------------------------------------------------------------------------*/

static Policy  cur_policy;
static size_t  threshold_T;       /* HTF only                               */
static size_t  heap_end;          /* bytes of `heap` managed by the fit list */
static size_t  rover;             /* Next Fit's resume point                 */

/* HTF bins: LIFO free lists of recycled small blocks.  The "next" pointer is
 * stored inside the free block's own payload - a free block's contents are
 * dead, so we may reuse them.  This is how real allocators do it too.      */
static size_t  bin_head[NUM_CLASSES];
#define BIN_NIL ((size_t)-1)

/* Buddy state */
static size_t  buddy_free[BUDDY_MAX_ORDER + 1];               /* list heads */
static unsigned char buddy_order_of[HEAP_SIZE >> BUDDY_MIN_ORDER];
static size_t  buddy_req_of[HEAP_SIZE >> BUDDY_MIN_ORDER];
static unsigned char buddy_is_free[HEAP_SIZE >> BUDDY_MIN_ORDER];

/* Measurement counters */
typedef struct {
    unsigned long long search_steps;   /* free-list nodes / orders examined  */
    unsigned long long worst_search;   /* worst single search                */
    unsigned long long n_alloc_calls;
    unsigned long long n_failures;
    size_t  live_bytes;                /* payload bytes currently handed out */
    size_t  requested_bytes;           /* what callers asked for, currently  */
    size_t  peak_live_requested;
    size_t  peak_internal_frag;
    size_t  bin_reserved;              /* HTF: bytes parked in bin lists     */
    size_t  peak_bin_reserved;
    double  peak_external_frag;
    unsigned long long validator_runs;
} Stats;

static Stats st;

/* ---------------------------------------------------------------------------
 * Initialisation
 * -------------------------------------------------------------------------*/

static void fit_heap_init(void)
{
    memset(heap, 0, HEAP_SIZE);
    heap_end = HEAP_SIZE;
    size_t payload = HEAP_SIZE - HDR - FTR;
    Header *h = hdr_at(0);
    h->size = payload; h->requested = 0; h->is_free = 1; h->bin_class = -1;
    ftr_at(0)->size = payload;
    rover = 0;
    for (int i = 0; i < NUM_CLASSES; i++) bin_head[i] = BIN_NIL;
}

static void buddy_heap_init(void)
{
    memset(heap, 0, HEAP_SIZE);
    memset(buddy_order_of, 0, sizeof buddy_order_of);
    memset(buddy_req_of, 0, sizeof buddy_req_of);
    memset(buddy_is_free, 0, sizeof buddy_is_free);
    for (int o = 0; o <= BUDDY_MAX_ORDER; o++) buddy_free[o] = BIN_NIL;
    buddy_free[BUDDY_MAX_ORDER] = 0;
    *(size_t *)(heap + 0) = BIN_NIL;          /* next pointer, in-block     */
    buddy_is_free[0] = 1;
    buddy_order_of[0] = BUDDY_MAX_ORDER;
}

static void heap_init(Policy p, size_t T)
{
    cur_policy = p;
    threshold_T = T;
    memset(&st, 0, sizeof st);
    if (p == BUDDY) buddy_heap_init(); else fit_heap_init();
}

/* ---------------------------------------------------------------------------
 * Fragmentation measurement
 *
 *  internal fragmentation = bytes inside live blocks the caller did not ask for
 *  external fragmentation = 1 - (largest free hole / total free bytes).
 *      0.0 means all free memory is in one usable piece; close to 1.0 means it
 *      is shattered into slivers.
 * -------------------------------------------------------------------------*/

static void sample_fragmentation(void)
{
    if (cur_policy == BUDDY) {
        size_t total_free = 0, largest = 0, internal = 0;
        for (int o = BUDDY_MIN_ORDER; o <= BUDDY_MAX_ORDER; o++) {
            for (size_t off = buddy_free[o]; off != BIN_NIL; off = *(size_t *)(heap + off)) {
                size_t sz = (size_t)1 << o;
                total_free += sz;
                if (sz > largest) largest = sz;
            }
        }
        for (size_t i = 0; i < (HEAP_SIZE >> BUDDY_MIN_ORDER); i++)
            if (buddy_order_of[i] && !buddy_is_free[i] && buddy_req_of[i])
                internal += ((size_t)1 << buddy_order_of[i]) - buddy_req_of[i];
        if (internal > st.peak_internal_frag) st.peak_internal_frag = internal;
        if (total_free) {
            double ef = 1.0 - (double)largest / (double)total_free;
            if (ef > st.peak_external_frag) st.peak_external_frag = ef;
        }
        return;
    }

    size_t total_free = 0, largest = 0, internal = 0, reserved = 0;
    for (size_t off = 0; off < heap_end; off = next_off(off)) {
        Header *h = hdr_at(off);
        if (h->bin_class >= 0) {
            if (h->is_free) reserved += h->size;        /* parked in a bin  */
            else            internal += h->size - h->requested;
        } else if (h->is_free) {
            total_free += h->size;
            if (h->size > largest) largest = h->size;
        } else {
            internal += h->size - h->requested;
        }
    }
    if (internal > st.peak_internal_frag) st.peak_internal_frag = internal;
    if (reserved > st.peak_bin_reserved)  st.peak_bin_reserved  = reserved;
    st.bin_reserved = reserved;
    if (total_free) {
        double ef = 1.0 - (double)largest / (double)total_free;
        if (ef > st.peak_external_frag) st.peak_external_frag = ef;
    }
}

/* ---------------------------------------------------------------------------
 * THE VALIDATOR
 *
 *  Walks the entire heap and asserts the invariants that must hold at all
 *  times.  Running this after every single allocation and free is what makes
 *  the results trustworthy - a silently corrupt heap would still print a
 *  pretty table.
 *
 *   V1  block sizes plus tag overhead sum exactly to the heap size
 *   V2  every footer agrees with its header (boundary tags intact)
 *   V3  no two ordinary free blocks are adjacent (coalescing really ran)
 *   V4  payloads are 8-byte aligned and never below the minimum
 *   V5  a live block never reports requested > size
 *   V6  a bin block's size equals its class size
 * -------------------------------------------------------------------------*/

static int validate_heap(const char *where)
{
    if (cur_policy == BUDDY) return 1;     /* buddy has its own structure    */

    size_t total = 0, off = 0;
    int prev_free_ordinary = 0;

    while (off < heap_end) {
        Header *h = hdr_at(off);

        if (h->size != ftr_at(off)->size) {               /* V2 */
            fprintf(stderr, "[VALIDATE] %s: footer mismatch at %zu\n", where, off);
            return 0;
        }
        if (h->size < MIN_PAYLOAD || (h->size % ALIGNMENT)) {   /* V4 */
            fprintf(stderr, "[VALIDATE] %s: bad payload size %zu at %zu\n",
                    where, h->size, off);
            return 0;
        }
        if (!h->is_free && h->requested > h->size) {      /* V5 */
            fprintf(stderr, "[VALIDATE] %s: requested>size at %zu\n", where, off);
            return 0;
        }
        if (h->bin_class >= 0) {                          /* V6 */
            size_t need = align_up(class_size[h->bin_class] < MIN_PAYLOAD
                                   ? MIN_PAYLOAD : class_size[h->bin_class]);
            /* It may be LARGER than the class size: split_block deliberately
             * refuses to split when the leftover would be too small to use. */
            if (h->size < need) {
                fprintf(stderr, "[VALIDATE] %s: bin block smaller than its "
                                "class at %zu\n", where, off);
                return 0;
            }
        }

        int ordinary_free = (h->is_free && h->bin_class < 0);
        if (ordinary_free && prev_free_ordinary) {        /* V3 */
            fprintf(stderr, "[VALIDATE] %s: uncoalesced free pair at %zu\n", where, off);
            return 0;
        }
        prev_free_ordinary = ordinary_free;

        total += block_bytes(off);
        off = next_off(off);
    }

    if (total != heap_end) {                              /* V1 */
        fprintf(stderr, "[VALIDATE] %s: sizes sum to %zu, expected %zu\n",
                where, total, heap_end);
        return 0;
    }
    st.validator_runs++;
    return 1;
}

/* ---------------------------------------------------------------------------
 * Block splitting
 *
 *  Only split if the leftover can itself hold a usable block.  Without this
 *  rule, Best Fit grinds the heap into slivers nobody can ever allocate - the
 *  classic failure mode we are trying to measure and avoid.
 * -------------------------------------------------------------------------*/

static void split_block(size_t off, size_t need)
{
    Header *h = hdr_at(off);
    size_t leftover = h->size - need;
    if (leftover < HDR + FTR + MIN_PAYLOAD) return;        /* keep it whole  */

    size_t rest_payload = leftover - HDR - FTR;

    h->size = need;
    ftr_at(off)->size = need;

    size_t roff = off + HDR + need + FTR;
    Header *r = hdr_at(roff);
    r->size = rest_payload; r->requested = 0; r->is_free = 1; r->bin_class = -1;
    ftr_at(roff)->size = rest_payload;
}

/* ---------------------------------------------------------------------------
 * The four classical fit policies
 * -------------------------------------------------------------------------*/

static size_t find_fit(size_t need, Policy p, unsigned long long *steps)
{
    size_t best = BIN_NIL, best_size = 0;
    *steps = 0;

    if (p == NEXT_FIT) {
        if (rover >= heap_end) rover = 0;      /* defensive: never dangle */
        /* two passes: from the rover to the end, then wrap to the rover */
        for (int pass = 0; pass < 2; pass++) {
            size_t start = (pass == 0) ? rover : 0;
            size_t stop  = (pass == 0) ? heap_end : rover;
            for (size_t off = start; off < stop; off = next_off(off)) {
                (*steps)++;
                Header *h = hdr_at(off);
                if (h->is_free && h->bin_class < 0 && h->size >= need) {
                    rover = off;
                    return off;
                }
            }
        }
        return BIN_NIL;
    }

    for (size_t off = 0; off < heap_end; off = next_off(off)) {
        (*steps)++;
        Header *h = hdr_at(off);
        if (!h->is_free || h->bin_class >= 0 || h->size < need) continue;

        if (p == FIRST_FIT) return off;
        if (p == BEST_FIT) {
            if (best == BIN_NIL || h->size < best_size) { best = off; best_size = h->size; }
            if (h->size == need) break;                   /* perfect fit    */
        }
        if (p == WORST_FIT) {
            if (best == BIN_NIL || h->size > best_size) { best = off; best_size = h->size; }
        }
    }
    return best;
}

/* ---------------------------------------------------------------------------
 * Coalescing (boundary tags make the backward merge possible)
 * -------------------------------------------------------------------------*/

static void coalesce(size_t off)
{
    Header *h = hdr_at(off);
    if (h->bin_class >= 0) return;          /* bin blocks stay in their bin  */

    /* forward */
    size_t n = next_off(off);
    if (n < heap_end) {
        Header *nh = hdr_at(n);
        if (nh->is_free && nh->bin_class < 0) {
            /* Next Fit's rover may be sitting on the block we are about to
             * swallow.  Leaving it there would leave it pointing into the
             * middle of a block, and the next search would read garbage
             * tags.  This was found by the validator, not by guesswork. */
            if (rover == n) rover = off;
            h->size += block_bytes(n);
            ftr_at(off)->size = h->size;
        }
    }
    /* backward: read the footer sitting immediately before our header */
    if (off >= HDR + FTR + MIN_PAYLOAD) {
        Footer *pf = (Footer *)(heap + off - FTR);
        size_t prev_payload = pf->size;
        size_t poff = off - FTR - prev_payload - HDR;
        if (poff < off) {
            Header *ph = hdr_at(poff);
            if (ph->size == prev_payload && ph->is_free && ph->bin_class < 0) {
                ph->size += block_bytes(off);
                ftr_at(poff)->size = ph->size;
                if (rover == off) rover = poff;
            }
        }
    }
}

/* ---------------------------------------------------------------------------
 * Buddy allocator
 * -------------------------------------------------------------------------*/

static int order_for(size_t size)
{
    int o = BUDDY_MIN_ORDER;
    while (((size_t)1 << o) < size && o < BUDDY_MAX_ORDER) o++;
    return o;
}
static size_t bidx(size_t off) { return off >> BUDDY_MIN_ORDER; }

static void buddy_list_push(int o, size_t off)
{
    *(size_t *)(heap + off) = buddy_free[o];
    buddy_free[o] = off;
    buddy_is_free[bidx(off)] = 1;
    buddy_order_of[bidx(off)] = (unsigned char)o;
}

static int buddy_list_remove(int o, size_t target)
{
    size_t *link = &buddy_free[o];
    while (*link != BIN_NIL) {
        if (*link == target) { *link = *(size_t *)(heap + target); return 1; }
        link = (size_t *)(heap + *link);
    }
    return 0;
}

static void *buddy_alloc(size_t size, unsigned long long *steps)
{
    *steps = 0;
    if (size == 0) return NULL;
    int want = order_for(size);
    int o = want;
    while (o <= BUDDY_MAX_ORDER && buddy_free[o] == BIN_NIL) { (*steps)++; o++; }
    if (o > BUDDY_MAX_ORDER) return NULL;

    size_t off = buddy_free[o];
    buddy_free[o] = *(size_t *)(heap + off);
    buddy_is_free[bidx(off)] = 0;

    while (o > want) {                       /* split down to the wanted order */
        o--;
        size_t buddy = off + ((size_t)1 << o);
        buddy_list_push(o, buddy);
        (*steps)++;
    }
    buddy_order_of[bidx(off)] = (unsigned char)want;
    buddy_req_of[bidx(off)]   = size;
    buddy_is_free[bidx(off)]  = 0;
    return heap + off;
}

static void buddy_free_block(void *p)
{
    size_t off = (size_t)((unsigned char *)p - heap);
    int o = buddy_order_of[bidx(off)];
    buddy_req_of[bidx(off)] = 0;

    while (o < BUDDY_MAX_ORDER) {
        size_t buddy = off ^ ((size_t)1 << o);
        if (!buddy_is_free[bidx(buddy)] || buddy_order_of[bidx(buddy)] != o) break;
        if (!buddy_list_remove(o, buddy)) break;
        buddy_is_free[bidx(buddy)] = 0;
        if (buddy < off) off = buddy;        /* merged block starts lower     */
        o++;
    }
    buddy_list_push(o, off);
}

/* ---------------------------------------------------------------------------
 *                        HYBRID THRESHOLD FIT  (HTF)
 *
 *  HTF_ALLOCATE(size):
 *    1  if size <= T:
 *    2       c <- smallest size class >= size
 *    3       if bin[c] non-empty: pop and return it            -> O(1)
 *    4       else carve a fresh block of class c (First Fit) and tag it
 *    5  else:
 *    6       pick the smallest free block >= size              -> Best Fit
 *    7       split it if the leftover is usable
 *    8  if nothing found: return NULL
 *
 *  HTF_FREE(block):
 *    1  if the block is bin-owned: push it back onto its bin, do NOT coalesce
 *    2  else mark free and merge with the free neighbours on both sides
 *
 *  THE IDEA: small allocations never touch the large free region, so the
 *  thousands of tiny short-lived objects can no longer shatter the big holes
 *  that large allocations depend on.  Recycled bin blocks are reused as-is,
 *  so a small alloc/free cycle costs no searching at all after warm-up.
 *
 *  THE PRICE, stated honestly: a block parked in a bin is reserved for that
 *  size class and is not returned to the general pool.  If the workload stops
 *  using that class, those bytes sit idle.  The program measures this as
 *  "bin reserved" so the trade-off appears in the results instead of being
 *  hidden.  This is why the threshold T matters: too high, and HTF reserves
 *  too much; too low, and small objects go back to fragmenting the heap.
 * -------------------------------------------------------------------------*/

static int class_for(size_t size)
{
    for (int i = 0; i < NUM_CLASSES; i++) if (class_size[i] >= size) return i;
    return -1;
}

static void *htf_alloc(size_t size, unsigned long long *steps)
{
    *steps = 0;
    if (size == 0) return NULL;

    if (size <= threshold_T) {
        int c = class_for(size);
        if (c >= 0) {
            if (bin_head[c] != BIN_NIL) {                 /* O(1) recycle    */
                size_t off = bin_head[c];
                bin_head[c] = *(size_t *)payload_at(off);
                Header *h = hdr_at(off);
                h->is_free = 0;
                h->requested = size;
                *steps = 1;
                return payload_at(off);
            }
            size_t need = align_up(class_size[c] < MIN_PAYLOAD ? MIN_PAYLOAD
                                                               : class_size[c]);
            size_t off = find_fit(need, FIRST_FIT, steps);
            if (off == BIN_NIL) return NULL;
            split_block(off, need);
            Header *h = hdr_at(off);
            h->is_free = 0; h->requested = size; h->bin_class = c;
            return payload_at(off);
        }
    }

    size_t need = align_up(size < MIN_PAYLOAD ? MIN_PAYLOAD : size);
    size_t off = find_fit(need, BEST_FIT, steps);
    if (off == BIN_NIL) return NULL;
    split_block(off, need);
    Header *h = hdr_at(off);
    h->is_free = 0; h->requested = size; h->bin_class = -1;
    return payload_at(off);
}

static void htf_free(void *p)
{
    size_t off = off_of_payload(p);
    Header *h = hdr_at(off);
    h->is_free = 1;
    h->requested = 0;
    if (h->bin_class >= 0) {
        *(size_t *)payload_at(off) = bin_head[h->bin_class];
        bin_head[h->bin_class] = off;
        return;
    }
    coalesce(off);
}

/* ---------------------------------------------------------------------------
 * Public allocator interface used by every workload
 * -------------------------------------------------------------------------*/

#define VALIDATE_EVERY_OP 1

static void *ma(size_t size)
{
    unsigned long long steps = 0;
    void *p = NULL;

    st.n_alloc_calls++;

    if (cur_policy == BUDDY)      p = buddy_alloc(size, &steps);
    else if (cur_policy == HTF)   p = htf_alloc(size, &steps);
    else {
        size_t need = align_up(size < MIN_PAYLOAD ? MIN_PAYLOAD : size);
        size_t off = find_fit(need, cur_policy, &steps);
        if (off != BIN_NIL) {
            split_block(off, need);
            Header *h = hdr_at(off);
            h->is_free = 0; h->requested = size; h->bin_class = -1;
            p = payload_at(off);
        }
    }

    st.search_steps += steps;
    if (steps > st.worst_search) st.worst_search = steps;
    if (!p) { st.n_failures++; return NULL; }

    st.requested_bytes += size;
    if (st.requested_bytes > st.peak_live_requested)
        st.peak_live_requested = st.requested_bytes;

    sample_fragmentation();
#if VALIDATE_EVERY_OP
    if (!validate_heap("after alloc")) exit(2);
#endif
    return p;
}

static void mf(void *p, size_t size)
{
    if (!p) return;
    if (cur_policy == BUDDY) buddy_free_block(p);
    else if (cur_policy == HTF) htf_free(p);
    else {
        size_t off = off_of_payload(p);
        hdr_at(off)->is_free = 1;
        hdr_at(off)->requested = 0;
        coalesce(off);
    }
    st.requested_bytes -= size;
    sample_fragmentation();
#if VALIDATE_EVERY_OP
    if (!validate_heap("after free")) exit(2);
#endif
}

/* ---------------------------------------------------------------------------
 * Deterministic RNG, so every run reproduces exactly (xorshift32)
 * -------------------------------------------------------------------------*/

static uint32_t rng_state = 2463534242u;
static void     rng_seed(uint32_t s) { rng_state = s ? s : 1u; }
static uint32_t rng_next(void)
{
    uint32_t x = rng_state;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    return rng_state = x;
}
static size_t rng_range(size_t lo, size_t hi) /* inclusive */
{
    return lo + (size_t)(rng_next() % (uint32_t)(hi - lo + 1));
}

/* ---------------------------------------------------------------------------
 * Synthetic workloads
 *
 *  Each one is a fixed sequence of allocate/free operations.  The SAME
 *  sequence (same seed) is replayed for every policy, so the comparison is
 *  apples to apples.
 * -------------------------------------------------------------------------*/

#define MAX_LIVE 20000

typedef struct { void *p; size_t size; } Slot;
static Slot live[MAX_LIVE];
static int  n_live;

static void live_reset(void) { n_live = 0; }

static void live_add(void *p, size_t size)
{
    if (!p || n_live >= MAX_LIVE) { if (p) mf(p, size); return; }
    live[n_live].p = p; live[n_live].size = size; n_live++;
}

static void live_free_index(int i)
{
    if (i < 0 || i >= n_live) return;
    mf(live[i].p, live[i].size);
    live[i] = live[n_live - 1];
    n_live--;
}

/* W1 - BURST: many small, short-lived objects.  The common case. */
static void workload_burst(void)
{
    live_reset(); rng_seed(11);
    for (int i = 0; i < 10000; i++) {
        size_t s = rng_range(8, 128);
        live_add(ma(s), s);
        if (n_live > 200) live_free_index((int)(rng_next() % (uint32_t)n_live));
    }
    while (n_live) live_free_index(n_live - 1);
}

/* W2 - ACCUMULATION: large, long-lived objects that are rarely released. */
static void workload_accumulation(void)
{
    live_reset(); rng_seed(22);
    for (int i = 0; i < 2000; i++) {
        size_t s = rng_range(1024, 8192);
        live_add(ma(s), s);
        if (n_live > 60) live_free_index((int)(rng_next() % (uint32_t)n_live));
    }
    while (n_live) live_free_index(n_live - 1);
}

/* W3 - BIMODAL: the realistic mix.  85% small and short-lived, 15% large and
 * long-lived.  This is the workload the threshold T exists for. */
static void workload_bimodal(void)
{
    live_reset(); rng_seed(33);
    for (int i = 0; i < 8000; i++) {
        size_t s = (rng_next() % 100u < 85u) ? rng_range(8, 128)
                                             : rng_range(1024, 16384);
        live_add(ma(s), s);
        if (n_live > 150) live_free_index((int)(rng_next() % (uint32_t)n_live));
    }
    while (n_live) live_free_index(n_live - 1);
}

/* W4 - ADVERSARIAL: built on purpose to shatter the heap.  Allocate an
 * alternating pattern of two sizes, free every other block, then ask for
 * something large.  Plenty of memory is free, but in slivers - so the request
 * fails.  Being able to construct the worst case proves you understand the
 * mechanism rather than merely coding it. */
static void workload_adversarial(void)
{
    live_reset(); rng_seed(44);
    const int N = 1200;
    for (int i = 0; i < N; i++) {
        size_t s = (i % 2 == 0) ? 64 : 512;
        live_add(ma(s), s);
    }
    for (int i = n_live - 1; i >= 0; i--)          /* free every other one   */
        if (i % 2 == 0) live_free_index(i);
    for (int i = 0; i < 40; i++) {                 /* now ask for big blocks */
        void *p = ma(32768);
        if (p) live_add(p, 32768);
    }
    while (n_live) live_free_index(n_live - 1);
}

/* ---------------------------------------------------------------------------
 * REAL PROGRAMS running on the allocator
 *
 *  These are not traces - they are ordinary C programs using real pointers,
 *  real structs and real strings, whose malloc/free go through the allocator
 *  under test.  If the allocator were wrong, these would produce wrong
 *  answers or crash.  They check correctness as well as performance.
 * -------------------------------------------------------------------------*/

/* R1 - linked list: build, delete every other node, then extend. */
typedef struct Node { int value; char *label; struct Node *next; } Node;

static long real_linked_list(void)
{
    Node *head = NULL;
    long checksum = 0;

    for (int i = 0; i < 3000; i++) {
        Node *n = (Node *)ma(sizeof(Node));
        if (!n) break;
        n->label = (char *)ma(24);
        if (n->label) snprintf(n->label, 24, "node-%d", i);
        n->value = i;
        n->next = head;
        head = n;
    }
    /* delete every other node */
    Node *prev = NULL, *cur = head;
    int k = 0;
    while (cur) {
        Node *nxt = cur->next;
        if (k % 2 == 0) {
            if (prev) prev->next = nxt; else head = nxt;
            if (cur->label) mf(cur->label, 24);
            mf(cur, sizeof(Node));
        } else prev = cur;
        cur = nxt; k++;
    }
    /* extend again, reusing the holes just created */
    for (int i = 0; i < 1500; i++) {
        Node *n = (Node *)ma(sizeof(Node));
        if (!n) break;
        n->label = (char *)ma(24);
        if (n->label) snprintf(n->label, 24, "extra-%d", i);
        n->value = 100000 + i;
        n->next = head;
        head = n;
    }
    for (Node *p = head; p; p = p->next) checksum += p->value + (long)strlen(p->label ? p->label : "");
    while (head) {
        Node *nxt = head->next;
        if (head->label) mf(head->label, 24);
        mf(head, sizeof(Node));
        head = nxt;
    }
    return checksum;
}

/* R2 - text processing: split a generated sentence stream into words, copy
 * each word onto the heap, keep a table, then release the duplicates. */
static long real_text_processing(void)
{
    static const char *words[] = {
        "operating", "system", "memory", "allocator", "fragmentation",
        "threshold", "coalescing", "boundary", "tag", "heap"
    };
    const int CAP = 4000;
    char  **table = (char **)ma(sizeof(char *) * CAP);
    size_t *sizes = (size_t *)ma(sizeof(size_t) * CAP);
    if (!table || !sizes) return 0;
    int n = 0; long checksum = 0;

    for (int i = 0; i < CAP; i++) {
        const char *w = words[rng_next() % 10u];
        size_t len = strlen(w) + 8;
        char *copy = (char *)ma(len);
        if (!copy) break;
        snprintf(copy, len, "%s%d", w, i % 97);
        table[n] = copy; sizes[n] = len; n++;
    }
    for (int i = 0; i < n; i++) checksum += (long)strlen(table[i]);
    for (int i = 0; i < n; i++) mf(table[i], sizes[i]);
    mf(sizes, sizeof(size_t) * CAP);
    mf(table, sizeof(char *) * CAP);
    return checksum;
}

/* R3 - matrix: row-pointer array plus one allocation per row, multiply, free. */
static long real_matrix(void)
{
    const int N = 80;
    double **A = (double **)ma(sizeof(double *) * N);
    double **B = (double **)ma(sizeof(double *) * N);
    double **C = (double **)ma(sizeof(double *) * N);
    if (!A || !B || !C) return 0;

    for (int i = 0; i < N; i++) {
        A[i] = (double *)ma(sizeof(double) * N);
        B[i] = (double *)ma(sizeof(double) * N);
        C[i] = (double *)ma(sizeof(double) * N);
        if (!A[i] || !B[i] || !C[i]) return 0;
        for (int j = 0; j < N; j++) { A[i][j] = i + j; B[i][j] = (i == j); C[i][j] = 0; }
    }
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            double s = 0;
            for (int k2 = 0; k2 < N; k2++) s += A[i][k2] * B[k2][j];
            C[i][j] = s;
        }
    long checksum = 0;
    for (int i = 0; i < N; i++) checksum += (long)C[i][i];   /* B is identity */

    for (int i = 0; i < N; i++) {
        mf(A[i], sizeof(double) * N);
        mf(B[i], sizeof(double) * N);
        mf(C[i], sizeof(double) * N);
    }
    mf(A, sizeof(double *) * N);
    mf(B, sizeof(double *) * N);
    mf(C, sizeof(double *) * N);
    return checksum;
}

/* ---------------------------------------------------------------------------
 * Reporting
 * -------------------------------------------------------------------------*/

typedef void (*WorkloadFn)(void);

typedef struct {
    const char *name;
    WorkloadFn  fn;
} Workload;

static const Workload workloads[] = {
    { "W1 Burst (small, short-lived)",   workload_burst        },
    { "W2 Accumulation (large, long)",   workload_accumulation },
    { "W3 Bimodal (realistic mix)",      workload_bimodal      },
    { "W4 Adversarial (worst case)",     workload_adversarial  },
};
#define N_WORKLOADS ((int)(sizeof workloads / sizeof workloads[0]))

static void print_header(const char *title)
{
    printf("\n%s\n", title);
    printf("%-14s %12s %12s %10s %10s %9s %9s\n",
           "Policy", "InternalFrag", "ExternalFrag", "AvgSteps",
           "MaxSteps", "Failures", "Reserved");
    printf("--------------------------------------------------------------"
           "-------------------\n");
}

static void print_row(const char *name)
{
    double avg = st.n_alloc_calls ? (double)st.search_steps / (double)st.n_alloc_calls : 0.0;
    printf("%-14s %11zu B %11.3f %10.2f %10llu %9llu %7zu B\n",
           name, st.peak_internal_frag, st.peak_external_frag, avg,
           st.worst_search, st.n_failures, st.peak_bin_reserved);
}

static FILE *csv = NULL;

static void csv_row(const char *workload, const char *policy, size_t T)
{
    if (!csv) return;
    double avg = st.n_alloc_calls ? (double)st.search_steps / (double)st.n_alloc_calls : 0.0;
    fprintf(csv, "%s,%s,%zu,%zu,%.6f,%.4f,%llu,%llu,%zu,%llu\n",
            workload, policy, T, st.peak_internal_frag, st.peak_external_frag,
            avg, st.worst_search, st.n_failures, st.peak_bin_reserved,
            st.validator_runs);
}

/* ---------------------------------------------------------------------------
 * main
 * -------------------------------------------------------------------------*/

int main(int argc, char **argv)
{
    int want_csv = (argc > 1 && strcmp(argv[1], "--csv") == 0);
    if (want_csv) {
        csv = fopen("results.csv", "w");
        if (csv) fprintf(csv, "workload,policy,T,internal_frag,external_frag,"
                              "avg_steps,max_steps,failures,bin_reserved,"
                              "validator_runs\n");
    }

    printf("=================================================================="
           "=================\n");
    printf(" Hybrid Threshold Fit (HTF) - comparative evaluation\n");
    printf(" Heap %u bytes | header %zu B + footer %zu B per block | "
           "alignment %u B\n", HEAP_SIZE, HDR, FTR, ALIGNMENT);
    printf(" Heap invariants are validated after EVERY allocation and free.\n");
    printf("=================================================================="
           "=================\n");

    /* ---- Part 1: baselines vs HTF on each workload ---------------------- */
    const size_t T_default = 256;

    for (int w = 0; w < N_WORKLOADS; w++) {
        print_header(workloads[w].name);
        for (int p = FIRST_FIT; p <= HTF; p++) {
            heap_init((Policy)p, T_default);
            workloads[w].fn();
            char label[32];
            if ((Policy)p == HTF) snprintf(label, sizeof label, "HTF T=%zu", T_default);
            else                  snprintf(label, sizeof label, "%s", policy_name[p]);
            print_row(label);
            csv_row(workloads[w].name, label, (Policy)p == HTF ? T_default : 0);
        }
    }

    /* ---- Part 2: threshold tuning --------------------------------------- */
    const size_t Ts[] = { 32, 64, 128, 256, 512, 1024 };
    for (int w = 0; w < N_WORKLOADS; w++) {
        char title[96];
        snprintf(title, sizeof title, "Threshold tuning - %s", workloads[w].name);
        print_header(title);
        for (size_t i = 0; i < sizeof Ts / sizeof Ts[0]; i++) {
            heap_init(HTF, Ts[i]);
            workloads[w].fn();
            char label[32];
            snprintf(label, sizeof label, "HTF T=%zu", Ts[i]);
            print_row(label);
            csv_row(workloads[w].name, label, Ts[i]);
        }
    }

    /* ---- Part 3: real programs ------------------------------------------ */
    printf("\nReal programs running on the allocator\n");
    printf("%-14s %-22s %12s %10s %10s %9s\n",
           "Policy", "Program", "InternalFrag", "AvgSteps", "MaxSteps", "Result");
    printf("--------------------------------------------------------------"
           "-------------------\n");

    for (int p = FIRST_FIT; p <= HTF; p++) {
        struct { const char *name; long (*fn)(void); } progs[] = {
            { "R1 linked list",     real_linked_list    },
            { "R2 text processing", real_text_processing },
            { "R3 matrix multiply", real_matrix          },
        };
        for (int i = 0; i < 3; i++) {
            heap_init((Policy)p, T_default);
            rng_seed(777);
            long result = progs[i].fn();
            double avg = st.n_alloc_calls
                       ? (double)st.search_steps / (double)st.n_alloc_calls : 0.0;
            printf("%-14s %-22s %11zu B %10.2f %10llu %9ld\n",
                   policy_name[p], progs[i].name, st.peak_internal_frag,
                   avg, st.worst_search, result);
        }
    }

    printf("\nAll heap invariants held throughout every run.\n");
    if (csv) { fclose(csv); printf("Wrote results.csv\n"); }
    return 0;
}
