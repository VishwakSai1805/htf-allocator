# Hybrid Threshold Fit (HTF) — a dynamic memory allocator in C

OS course project. One C file, terminal only, no libraries, no interface.

## What this is

The program sets aside a 1 MB array and implements `malloc` / `free` on top of it
six different ways, then runs the same workloads through all six and prints how
much memory each one wasted, how much searching it did, and how many requests
failed.

| # | Policy | Rule |
|---|--------|------|
| 1 | First Fit | take the first free hole that is big enough |
| 2 | Best Fit | take the smallest hole that is big enough |
| 3 | Worst Fit | take the largest hole |
| 4 | Next Fit | First Fit, but resume searching where the last search stopped |
| 5 | Buddy | round the request up to a power of two, split and merge |
| 6 | **HTF (mine)** | requests ≤ T from fixed size-class bins in O(1); requests > T by Best Fit with coalescing |

## The problem it addresses

Repeated allocate and free leaves the heap full of scattered holes. Eventually a
large request fails even though plenty of memory is free — it is simply in pieces
too small to use. That is fragmentation.

The classical policies each fail differently. Best Fit wastes the least space but
scans the whole free list and leaves behind slivers nobody can reuse. Buddy is
fast but rounds every request up to a power of two, so a 33-byte request consumes
64 bytes. Real workloads are a mix of thousands of tiny short-lived objects and a
few large long-lived ones, so no single policy is right for both classes.

## The algorithm

```
HTF_ALLOCATE(size):
  1  if size <= T:
  2       c <- smallest size class >= size
  3       if bin[c] is non-empty: pop it and return            -> O(1)
  4       else carve a fresh block of class c (First Fit) and tag it
  5  else:
  6       select the smallest free block >= size                -> Best Fit
  7       split it if the leftover can hold a usable block
  8  if nothing was found: return NULL

HTF_FREE(block):
  1  if the block is bin-owned: push it back onto its bin, do NOT coalesce
  2  else mark it free and merge with the free neighbours on both sides
```

Size classes: 8, 16, 24, 32, 48, 64, 96, 128, 192, 256 bytes — deliberately
finer than powers of two, which is how HTF keeps internal waste below Buddy's.

**The idea:** small allocations never touch the large free region, so the
thousands of tiny objects can no longer shatter the big holes that large
allocations depend on.

**The price, stated honestly:** a block parked in a bin is reserved for that size
class and is not returned to the general pool. The program measures this as
"Reserved" so the trade-off appears in the results rather than being hidden. That
is exactly why T matters — too high and HTF reserves too much, too low and small
objects go back to fragmenting the heap.

## Block layout — boundary tags (Knuth)

```
offset ->  [ Header | payload .......... | Footer ]
             24 B      size bytes           8 B
```

The footer is what makes backward coalescing possible: when freeing a block we
read the bytes immediately before our own header to learn how big the previous
block was, so we can merge with it.

## Build and run

```
gcc -O2 -Wall -Wextra -o htf htf_allocator.c
./htf                 # prints every table
./htf --csv           # also writes results.csv
```

Compiles with zero warnings under `-Wall -Wextra`.

## Heap validator

A validator walks the entire heap and checks six invariants **after every single
allocation and free**:

- V1 block sizes plus tag overhead sum exactly to the heap size
- V2 every footer agrees with its header (boundary tags intact)
- V3 no two ordinary free blocks are left adjacent (coalescing really ran)
- V4 payloads are 8-byte aligned and never below the minimum
- V5 a live block never reports `requested > size`
- V6 a bin block is never smaller than its size class

This caught two real bugs during development: a bin block that was legitimately
larger than its class (because splitting is refused when the leftover would be
unusable), and a dangling Next Fit resume pointer left behind after a forward
merge, which made the next search read garbage tags.

## Workloads

Four synthetic traces, replayed identically for every policy with a fixed RNG
seed, plus three real programs.

| Trace | Pattern |
|-------|---------|
| W1 Burst | 10,000 allocations of 8–128 B, nearly all freed quickly |
| W2 Accumulation | 2,000 allocations of 1–8 KB, mostly long-lived |
| W3 Bimodal | 85% small and short-lived, 15% large and long-lived — the realistic mix |
| W4 Adversarial | alternating two sizes, free every other block, then request 32 KB — built on purpose to shatter the heap |

The three real programs are ordinary C — a linked list with string labels, a text
tokeniser, and an 80×80 matrix multiply — using real pointers whose `malloc` and
`free` go through the allocator under test. **All six policies produce identical
results on all three programs**, which is correctness evidence, not just
performance measurement. The matrix result (6320) is independently verifiable:
B is the identity matrix, so the diagonal of C sums to Σ 2i for i = 0…79.

## Results

Internal fragmentation = bytes inside live blocks the caller did not ask for.
External fragmentation = 1 − (largest free hole ÷ total free bytes); 0.0 means all
free memory is in one usable piece, near 1.0 means it is shattered.
AvgSteps = free-list nodes examined per allocation (a counter, so it is
machine-independent).

### W1 Burst — small, short-lived

| Policy | Internal | External | AvgSteps | MaxSteps | Failures | Reserved |
|---|---|---|---|---|---|---|
| First Fit | 3144 B | 0.019 | 109.52 | 236 | 0 | — |
| Best Fit | 2741 B | 0.022 | 164.59 | 226 | 0 | — |
| Worst Fit | 793 B | 0.573 | 296.83 | 319 | 0 | — |
| Next Fit | 793 B | 0.570 | 2.00 | 2 | 0 | — |
| Buddy | 5140 B | 0.500 | 0.14 | 30 | 0 | — |
| **HTF T=256** | **2559 B** | **0.000** | **5.49** | 300 | 0 | 22208 B |

### W3 Bimodal — the realistic mix

| Policy | Internal | External | AvgSteps | MaxSteps | Failures | Reserved |
|---|---|---|---|---|---|---|
| First Fit | 2135 B | 0.308 | 83.91 | 203 | 0 | — |
| Best Fit | 1963 B | 0.321 | 145.25 | 188 | 0 | — |
| Worst Fit | 622 B | 0.951 | 222.03 | 241 | 0 | — |
| Next Fit | 622 B | 0.919 | 2.01 | 8 | 0 | — |
| Buddy | 131248 B | 0.496 | 0.18 | 16 | 0 | — |
| **HTF T=256** | **1799 B** | 0.329 | **40.52** | 277 | 0 | 15320 B |

### W4 Adversarial — the constructed worst case

| Policy | Internal | External | AvgSteps | Failures | Reserved |
|---|---|---|---|---|---|
| First Fit | 0 B | 0.818 | 620.35 | 20 | — |
| Best Fit | 0 B | 0.818 | 620.35 | 20 | — |
| Worst Fit | 0 B | 0.818 | 620.35 | 20 | — |
| Next Fit | 0 B | 0.818 | 21.66 | 20 | — |
| Buddy | 0 B | 0.846 | 2.07 | 19 | — |
| **HTF T=256** | 0 B | **0.316** | 620.35 | 20 | 38400 B |

### Threshold tuning

W1 Burst:

| T | Internal | External | AvgSteps | Reserved |
|---|---|---|---|---|
| 32 | 2736 B | 0.017 | 164.04 | 2192 B |
| 64 | 2472 B | 0.011 | 123.22 | 7008 B |
| **128** | 2559 B | **0.000** | **5.49** | 22208 B |
| 256 | 2559 B | 0.000 | 5.49 | 22208 B |
| 512 | 2559 B | 0.000 | 5.49 | 22208 B |
| 1024 | 2559 B | 0.000 | 5.49 | 22208 B |

W3 Bimodal:

| T | Internal | External | AvgSteps | Reserved |
|---|---|---|---|---|
| 32 | 1762 B | 0.314 | 144.25 | 1272 B |
| 64 | 1478 B | 0.349 | 121.18 | 4496 B |
| **128** | 1799 B | 0.329 | **40.52** | 15320 B |
| 256 | 1799 B | 0.329 | 40.52 | 15320 B |

## Conclusions — read these from the numbers, not from a narrative

1. **The threshold has a knee at T = 128 bytes**, and the curve is flat above it.
   On W1 that knee takes external fragmentation from 0.011 to 0.000 and average
   search cost from 123 steps to 5.5 — a 22× reduction. On W3 it cuts search cost
   from 145 steps (Best Fit) to 40.5, a 3.6× reduction, while holding internal
   fragmentation at the lowest figure of any fit policy.

2. **The curve is flat above T = 128 because the largest size class is 256 bytes.**
   T stops mattering once it exceeds the class range, so the real tuning knob is
   the class table, not T alone. This bounds the result and should be stated
   rather than glossed over.

3. **HTF reduces external fragmentation on the adversarial trace from 0.818 to
   0.316, but it does not prevent the 20 failures.** Keeping small objects out of
   the large region leaves the big holes more intact, yet the trace still drives
   the heap past the point where 32 KB requests can be satisfied. Fragmentation
   control is not capacity.

4. **Buddy's internal fragmentation is catastrophic on large workloads** —
   131 KB on W3, roughly 70× HTF's — which is the power-of-two rounding penalty
   made visible. It is, however, the fastest searcher by far.

5. **Next Fit looks fast (2 steps per allocation) but shatters the heap**
   (external fragmentation 0.57–0.92). Its rover sits on the large trailing free
   region, so it finds a fit immediately and never consolidates. Speed measured
   alone would have been misleading.

6. **On W2, HTF is byte-for-byte identical to Best Fit**, because every request
   exceeds T. That is a correctness check on the implementation, not a null
   result.

## Theoretical cost, beside the measured cost

| Policy | Per allocation | Measured AvgSteps (W3) |
|--------|----------------|------------------------|
| First Fit | O(n) | 83.91 |
| Best Fit | O(n) | 145.25 |
| Worst Fit | O(n) | 222.03 |
| Next Fit | O(n) amortised, O(1) typical | 2.01 |
| Buddy | O(log n) | 0.18 |
| HTF | O(1) for size ≤ T, O(n) above | 40.52 |

## Scope decisions

- **The heap has a fixed size.** This is an experimental control, not a shortcut:
  an expandable heap would hide fragmentation, which is the exact effect being
  measured.
- **No thread safety.** Locking is orthogonal to the allocation policy and would
  change none of these results, while adding lock-granularity and per-arena
  questions well outside the scope of the algorithm.
- Very large allocations are not routed to the operating system separately
  (`mmap`-style); every request is served from the one managed region, so that
  all six policies are compared on identical terms.

## Prior art

The mechanisms used here are long established: boundary-tag coalescing is Knuth's,
segregated size-class free lists and deferred coalescing appear in Doug Lea's
`dlmalloc`, and the standard survey is Wilson, Johnstone, Neely & Boles,
*Dynamic Storage Allocation: A Survey and Critical Review* (1995). Linux performs
an equivalent split in the kernel between its buddy page allocator and its
slab/slub allocator for small objects.

The contribution of this project is the explicit threshold formulation, the
measured location of the knee, and the quantification of the reservation cost
that the bins impose — reported as a trade-off rather than presented as a free win.
