# OpenMP loop patterns on AMD and NVIDIA

This subdirectory contains tests for OpenMP + do concurrent loop patterns for both AMD and NVIDIA.
The reason for this comparison is that we've ported MOM6 code using OpenMP + do concurrent for
NVIDIA only, and I would like to get a sense of how well our code would translate to AMD. Hence,
the tests are heavily reduced versions of some common loop patterns used so far.

Setups:

* AMD: MI250x single GCD using ROCm 7.2.3.
  * Run inside Singularity with ROCm 7.2.3 in the container.
  * Host environment: ROCm 6.4.1 module available, Linux kernel 6.4.
  * Code compiled with `amdflang -fopenmp --offload-arch=gfx90a -Mnofma -O3 -ffp-contract=off -fdo-concurrent-to-openmp=device`
  * 28 TFLOPS, 1.6TB/s
* NVIDIA: V100 SXM using NVHPC 25.9.
  * Code compiled with `nvfortran -mp=gpu -gpu=mem:separate -Mnofma -Mnovect -O4 -Minline=flux_elem -stdpar=gpu`
  * 7.5 TFLOPS, 0.9TB/s
* AMD: Radeon 780M integrated GPU using ROCm 7.2.3. Used only for the correctness notes
  near the end of this file, not for any of the timing tables.
  * Bare metal, no container. `gfx1103`, RDNA3, 12 CUs, no XNACK.
  * Code compiled with `amdflang -fopenmp --offload-arch=gfx1103 -O3 -ffp-contract=off -fdo-concurrent-to-openmp=device`
  * Note the missing `-Mnofma`: that is an NVHPC/PGI-lineage flag and `amdflang` rejects it
    outright. `-ffp-contract=off` is LLVM's way of forbidding the same FMA fusion.

The `-ffp-contract=off` and `-Mnofma -Mnovect` flags are to ensure that results are bitwise
identical with CPU. Interestingly, excluding these flags *slows down* the NVIDIA GPU code.
`-Minline=flux_elem` is included in the NVIDIA compilation as `amdflang` automatically inlines at `-O3`
(can be checked with `-Rpass=inline`). Furthermore, there is a compiler bug affecting some of these
tests which require inlining functions to get correct results.

It's important to note that MI250x is a few years newer than the V100, and a more appropriate
comparison would be with the A100s. But the latter were harder to get a hold of at the time.

All table entries are timings in milliseconds, reported as the best of 5 runs. The tables show
all tested grid sizes; the largest is 512x512x100. The corresponding source files are:
pattern 1 is `repro_tests.F90`, pattern 2 is `av_rem_test.F90`, and pattern 3 is
`col_norm_test.F90`.

Rows labelled "with thread spec" explicitly set the team/thread geometry. In these tests that
means using the target-specific constants in `omp_macros.inc`: `TEAM_SIZE=256` and
`WAVEFRONT=64` on AMD, and `TEAM_SIZE=128` and `WAVEFRONT=32` on NVIDIA. The jki version
uses `thread_limit(round_up(nx, WAVEFRONT))`; fused ij versions use
`num_teams(ceil(nx * ny / TEAM_SIZE))`. Here `round_up(nx, WAVEFRONT)` means the next
multiple of `WAVEFRONT` greater than or equal to `nx`.

Note that wherever the `loop` construct is used in the below examples, the AMD code actually uses
the older `teams distribute parallel do`. On the other hand, NVIDIA uses `loop`. This is controlled
by `omp_macros.inc`: defining `LOOP` maps the macros to `loop`, while leaving it undefined maps them
to `teams distribute parallel do`. There are two
reasons for this difference:

1. With `amdflang`, `loop` didn't always give the right answers or they were slow.
2. NVIDIA performs significantly faster with `loop` than the older construct.

## Loop pattern 1a: Outer target teams with parallel inner ij loops + serial loops

This is probably the most complex loop pattern used so far. The `target teams` region is started
and then an initialisation ij loop is run. Following is a newton-raphson iterative loop that ought
to be serialised for every thread. Inside the newton-raphson loop is more parallel ij loops, with
one of them being inside a serialised k loop. The k loop needs to be serialised because there are
sum reductions and bitwise reproducibility is a goal.

The general code looks like:

```fortran
!$omp target teams

!$omp loop collapse(2) ! or distribute parallel do
do j=... ; do i=...

! newton-raphson loop
do itt=1,20
  !$omp loop collapse(2)
  do j=... ; do i=...
  do k=...
    !$omp loop collapse(2)
    do j=... ; do i=...
  !$omp loop collapse(2)
  do j=... ; do i=...
```

|                   | 32x32x100 | 64x64x100 | 128x128x100 | 256x256x100 | 512x512x100 |
| :---              | ---:      | ---:      | ---:        | ---:        | ---:        |
| MI250x ROCm 7.2.3 | 1.948     | 2.154     | 2.173       | 3.489       | 13.620      |
| V100 NVHPC 25.9   | 1.168     | 1.454     | 1.658       | 4.089       | 12.060      |

Despite being older and weaker, the V100 beats the MI250x in most problem sizes; the MI250x is
faster for the 256x256x100 case.

In this scenario, I couldn't get the `loop` construct working on the inner loops in the AMD setup.

Worth noting that when using ROCm 6.4.1, the AMD setup was about 10-15% slower than the above reported
times.

## Loop pattern 1b: Outer ij loop with inner serial loops

The function of this version of the code is the same as above, except structured in a more GPU-friendly
way. Instead of the ij loops being inside, the ij loop is outside with inner loops only being the serial
loops. This is beneficial because it is much clearer to the compiler that each GPU thread works
independently of the others, whereas in version 1a, the compiler had to infer. 

The code looks like:

```fortran
!$omp target teams loop collapse(2) ! or teams distribute parallel do
do j=... ; do i=...
  do itt=1,20
    do k=...
```

Which is clearly simpler.

|                   | 32x32x100 | 64x64x100 | 128x128x100 | 256x256x100 | 512x512x100 |
| :---              | ---:      | ---:      | ---:        | ---:        | ---:        |
| MI250x ROCm 7.2.3 | 0.548     | 0.697     | 0.722       | 1.496       | 6.671       |
| V100 NVHPC 25.9   | 0.374     | 0.667     | 0.859       | 2.234       | 8.881       |

And in both cases, the compiler does a better job of optimising. Furthermore, the AMD setup is now
doing better than the tested NVIDIA setup by a meaningful margin at larger problem sizes.

These timings didn't change much with ROCm version used.

## Loop pattern 2: k-reduction

This loop pattern is related to, but simpler than the first loop pattern. This simply reduces a 3d
array by the 3rd index i.e. ijk -> ij. Importantly, because bitwise reproducibility is an aim, the
reduction of the 3rd index must be serialised. In the ported MOM6 code, this is a jki loop to
preserve CPU performance.

```fortran
!$omp target teams loop ! or target teams distribute
do j=...
  do k=... ! must be serialised
    !$omp loop ! or parallel do
    do i=...
```

|                                    | 32x32x100 | 64x64x100 | 128x128x100 | 256x256x100 | 512x512x100 |
| :---                               | ---:      | ---:      | ---:        | ---:        | ---:        |
| MI250x ROCm 6.4.1                  | 0.956     | 1.823     | 3.713       | 7.440       | 16.536      |
| MI250x ROCm 7.2.3                  | 0.109     | 0.111     | 0.127       | 0.193       | 1.217       |
| MI250x ROCm 7.2.3 with thread spec | 0.103     | 0.106     | 0.123       | 0.300       | 0.759       |
| V100 NVHPC 25.9                    | 0.038     | 0.058     | 0.071       | 0.169       | 0.681       |
| V100 NVHPC 25.9 with thread spec   | 0.038     | 0.055     | 0.069       | 0.142       | 0.562       |

Notably, **when using ROCm 6.4.1, these timings were an order of magnitude worse than 7.2.3**,
which highlights how the AMD compilers are improving. 

Even with the newer ROCm, the AMD setup is slower than NVIDIA by approx 2x - except for the
256x256x100 test. This might be because the block size used in AMD is 256, which matches up 
nicely with that problem size.

I did notice that the 512x512x100 was reduced to ~0.8ms if the number of threads was set to match
the `nx` dimension. In contrast, the 256x256x100 time doubled, and the smaller sizes
improved a negligible amount. Now, given that problem sizes are unlikely to line up with
the default block size (256), it's probably beneficial to set threads explicitly when using
this jki loop pattern.

When using `thread_limit` in the NVIDIA setup, which uses the `loop` clause, the compiler
feedback indicated i-loops after the first one were being serialised:

```
# without thread_limit
         57, Generating "nvkernel_av_rem_mod_run_av_rem_omp__F1L57_2" GPU kernel
             Generating NVIDIA GPU code
           59, Loop parallelized across teams ! blockidx%x
           61, Loop parallelized across threads(128) ! threadidx%x
           64, Loop run sequentially 
           66, Loop parallelized across threads(128) ! threadidx%x

# with thread_limit
         57, Generating "nvkernel_av_rem_mod_run_av_rem_omp__F1L57_2" GPU kernel
             Generating NVIDIA GPU code
           59, Loop parallelized across teams ! blockidx%x
           61, Loop parallelized across threads(nthreads) ! threadidx%x
           64, Loop run sequentially 
           66, Loop run sequentially
```

But as the table shows, there was nevertheless a speedup. Perhaps a bug in the compiler info?


If turned into a jik loop:

```fortran
!$omp target teams loop collapse(2)
do j=...
  do i=...
    do k=...
```

then we get:

|                   | 32x32x100 | 64x64x100 | 128x128x100 | 256x256x100 | 512x512x100 |
| :---              | ---:      | ---:      | ---:        | ---:        | ---:        |
| MI250x ROCm 6.4.1 | 0.057     | 0.057     | 0.059       | 0.108       | 0.382       |
| MI250x ROCm 7.2.3 | 0.051     | 0.052     | 0.054       | 0.093       | 0.349       |
| V100 NVHPC 25.9   | 0.014     | 0.016     | 0.041       | 0.140       | 0.541       |

The timings look like the difference between loop pattern 1a and 1b. There was a small
improvement when using the newer ROCm.

### Do concurrent

I also tested do concurrent for this loop. It's worth noting that ROCm 6.4.1 doesn't
support `do concurrent` offload, but ROCm 7+ does.

For the jki loop:

|                   | 32x32x100 | 64x64x100 | 128x128x100 | 256x256x100 | 512x512x100 |
| :---              | ---:      | ---:      | ---:        | ---:        | ---:        |
| MI250x ROCm 7.2.3 | 1.433     | 2.770     | 5.247       | 10.704      | 23.818      |
| V100 NVHPC 25.9   | 0.019     | 0.021     | 0.045       | 0.140       | 0.525       |

Which is clearly awful on the AMD setup. By comparison, the NVIDIA setup is close to
the OpenMP times. When inspecting the kernel launches on the AMD setup with the
environment variable `LIBOMPTARGET_KERNEL_TRACE=1`, the output shows that the kernel is
being launched with 16 blocks, each with 32 threads for all the problem sizes, when
we would hope for `nj` blocks and 256 threads.

Interestingly, if we rerun the OpenMP jki loop version, except remove the inner
`parallel do`s and have a plain `target teams distribute parallel do` wrap the outer
j loop, we get the same times and launch configurations here. So we can infer that
`do concurrent(j=...)` is equivalent to
```
!$omp target teams distribute parallel do
do j=...
```
which is undesirable, as the compiler isn't parallelizing over i as well.

If we run the jik do concurrent version:

|                   | 32x32x100 | 64x64x100 | 128x128x100 | 256x256x100 | 512x512x100 |
| :---              | ---:      | ---:      | ---:        | ---:        | ---:        |
| MI250x ROCm 7.2.3 | 0.878     | 1.731     | 6.411       | 65.662      | 151.468     |
| V100 NVHPC 25.9   | 0.016     | 0.018     | 0.042       | 0.139       | 0.524       |

The NVIDIA timings are similar to previous, but the AMD setup seems to perform much
worse. Inspecting the launch configuration, we see that we're getting 440 blocks and
256 threads - which is what we'd hope, but this means the slowdown is for other reasons.

## Loop pattern 3: repeated simple loops

By simple loops, I mean 2d or 3d loops where each iteration is independent and all
the iterations can be simply parallelised. This test is useful to examine if `amdflang`'s
async OpenMP kernel launches proves useful or not.

```fortran
!$omp target loop ! or distribute parallel do
do j=... ; do i=...

do k=...
  !$omp target loop
  do j=... ; do i=...

!$omp target loop
do j=... ; do i=...

!$omp target loop
do k=... ; do j=... ; do i=...
```

The 2nd loop where the k is outside is intentional, to exaggerate the effect of kernel
launch overheads.

|                   | 32x32x100 | 64x64x100 | 128x128x100 | 256x256x100 | 512x512x100 |
| :---              | ---:      | ---:      | ---:        | ---:        | ---:        |
| MI250x ROCm 7.2.3 | 1.530     | 1.543     | 1.403       | 1.705       | 2.229       |
| V100 NVHPC 25.9   | 0.863     | 0.849     | 0.883       | 1.066       | 1.931       |

This appears to highlight that launch overheads are much higher on the AMD setup. And
despite the asynchronous kernel launches, the async launches cannot hide the launch
overhead for many small kernels.

## NVHPC 26.3 directive-form resource usage and timings

I also built `repro_tests.F90` three ways with NVHPC 26.3 to compare the generated
kernel resource usage for the same macro-based source:

```sh
FC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3/compilers/bin/nvfortran
CUOBJDUMP=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3/cuda/13.1/bin/cuobjdump
BASE="-mp=gpu -acc=gpu -gpu=mem:separate -O4 -stdpar=gpu -Minline=name:flux_elem -Mnovect -Mnofma"

make -B repro_tests FC="$FC" FFLAGS="$BASE"        # OpenMP long form
make -B repro_tests FC="$FC" FFLAGS="-DLOOP $BASE" # OpenMP loop
make -B repro_tests FC="$FC" FFLAGS="-DACC $BASE"  # OpenACC
```

The results below were collected with `run_nvhpc_directive_comparison.sh` in
`nvhpc_directive_results/20260522_125745` on an NVIDIA GeForce RTX 4060 Laptop GPU
(`sm_89`, driver 580.142). All three binaries exited successfully.

The compiler emitted cubins for several SM targets; the table below uses the `sm_89`
entries from `cuobjdump -res-usage`. Entries are `REG/STACK/SHARED`. For the long-form
OpenMP case, `cuobjdump` also reports outlined helper functions for the inner loops;
the table uses the enclosing launched kernel entry.

| Kernel / region | OpenMP long form | OpenMP loop | OpenACC |
| :--- | ---: | ---: | ---: |
| Continuity, line 72 | 148 / 920 / 0 | 46 / 0 / 0 | 49 / 0 / 0 |
| Flux adjust init, line 162 | 40 / 0 / 0 | 40 / 0 / 0 | 38 / 0 / 0 |
| Flux adjust separate loops, line 174 | 148 / 1008 / 0 | 78 / 0 / 24 | 88 / 0 / 24 |
| Tiled flux adjust, line 311 | 148 / 1016 / 0 | 52 / 0 / 48 | 50 / 0 / 48 |
| IJ-outer scalar-private, line 578 | 82 / 0 / 0 | 90 / 0 / 0 | 90 / 0 / 0 |
| Fused IJ, line 693 | 148 / 1608 / 0 | 80 / 0 / 24 | 78 / 0 / 24 |
| J-outer / I-parallel, line 821 | 148 / 1080 / 0 | 95 / 0 / 24 | 93 / 0 / 24 |

The main pattern is that NVHPC's long-form OpenMP lowering is much heavier here:
several kernels hit around 148 registers and also require about 1-1.6 KB of stack
per thread. The `loop` and OpenACC forms are close to one another, use no stack in
these kernels, and differ by only a few registers for most regions.

For the largest case, 512x512x100, the best-of-5 timings were:

| Kernel / region | OpenMP long form (ms) | OpenMP loop (ms) | OpenACC (ms) |
| :--- | ---: | ---: | ---: |
| Continuity | 16.804 | 6.239 | 6.246 |
| Flux adjust separate loops | 121.425 | 32.596 | 30.156 |
| IJ-outer scalar-private | 27.070 | 28.968 | 28.972 |
| Fused IJ | 40.950 | 29.999 | 29.029 |
| J-outer / I-parallel | 63.488 | 32.921 | 31.444 |
| Tiled flux adjust | 107.388 | 33.014 | 31.305 |

The timing story matches the resource-usage story. `loop` is 1.4-3.7x faster than
long-form OpenMP for every 512x512x100 kernel except the scalar-private IJ-outer
variant, where long-form OpenMP is about 7% faster. OpenACC is effectively tied with
OpenMP `loop` for continuity and the IJ-outer variant, and is otherwise a few percent
faster than `loop` in this run.

## Correctness traps found while porting MOM6 to AMD

The sections above are all about speed. These three are about getting the right answer at
all, and were found while porting `MOM_vert_friction.F90` to `amdflang` on the Radeon 780M
setup listed at the top. Each one cost days, and each one is silent: nothing in the build
log, the run log, or the numbers themselves tells you anything is wrong.

The corresponding reproducers are `dc_exit_test.F90`, `map_components_test.F90` and
`private_stack_test.F90`, all built by `make`. They are deliberately tiny -- the point of
each is to be the smallest thing that still misbehaves.

### `do concurrent` silently stays on the host if an inner loop has an `exit`

`-fdo-concurrent-to-openmp=device` will not lower a `do concurrent` loop whose body
contains a nested ordinary loop with an `exit` statement. It compiles as host-only serial
code, with no error, no warning, and no remark even under `-Rpass`.

This is far nastier than it sounds, because the host code is by construction correct. The
answers match the CPU baseline bitwise, so **every correctness check passes** and the loop
looks ported. It is a false positive: what you have verified is that the host fallback
reproduces the host answer.

```fortran
! stays on the host
do concurrent (j=..., i=...)
  do k = 1, nz
    acc = acc + h(i,j,k) * k
    depth = depth + h(i,j,k)
    if (depth >= limit) exit
  enddo
enddo

! dispatches: same arithmetic, no branch leaves the loop early
do concurrent (j=..., i=...)
  do k = 1, nz
    if (depth < limit) then
      acc = acc + h(i,j,k) * k
      depth = depth + h(i,j,k)
    endif
  enddo
enddo
```

`./dc_exit` reports where each variant actually ran:

```
  PASS  do concurrent, inner exit
        ran on: HOST (no kernel was dispatched)
  PASS  do concurrent, exit rewritten as if
        ran on: DEVICE
  PASS  !$omp target, inner exit
        ran on: DEVICE
```

Three PASSes, one of which never touched the GPU. Note the third line: the *same* `exit`,
left completely untouched under a hand-written `!$omp target teams distribute parallel do`,
dispatches fine. So this is specific to the do-concurrent lowering pass, not to OpenMP or
to the hardware. Either rewrite the `exit` away or use an explicit directive.

The reproducer detects this from inside the loop with `omp_is_initial_device()`, which has
to be redeclared as a `pure` C binding to get past `do concurrent`'s purity rule. Outside a
test program the practical check is `LIBOMPTARGET_INFO=16` (or `LIBOMPTARGET_KERNEL_TRACE=1`,
as used further up this file): look for a `Launching kernel` line whose symbol name carries
the routine and line number of the loop you just ported, and treat its absence as a failure
regardless of what the numbers say.

### Derived-type array components are not mapped for you

On NVIDIA this whole class of problem is invisible, because `-gpu=mem:managed` /
`mem:unified` lets a kernel dereference anything. `gfx1103` has no XNACK and consumer RDNA
parts generally do not, so there is no safety net: every component a kernel touches has to
be named in a map clause.

Two separate things are not automatic, and `./map_components` demonstrates both:

| case | map clause | result |
| :--- | :--- | :--- |
| 1 | `map(to: st)`, components read in a `declare target` helper | **fault** |
| 2 | `map(to: st, st%alloc, st%ptr)`, same helper | works |
| 3 | no map clause at all, components read inline | **fault** |
| 4 | `target enter data` / `exit data` writeback | works |

Case 3 is what you get by porting a managed-memory nvfortran kernel and simply deleting the
map clauses. A derived-type dummy argument *is* mapped implicitly, but that map carries only
the struct; the arrays its descriptors point at are not attached, so the first dereference
faults at a host address.

Case 1 is the one that is easy to miss even when you are being careful, because nothing at
the map site names the offending component:

```fortran
!$omp target teams distribute parallel do collapse(2) map(to: st)
do j = ... ; do i = ...
  call read_components(st, i, j, out(i,j))   ! reads st%alloc, st%ptr
enddo ; enddo
```

This was the last bug in the MOM6 port. `vertvisc_coef`'s offloaded regions mapped `visc`
and then read `visc%Kv_shear` inside `find_coupling_coef_k`, several call levels down. It
faulted or hung outright, depending on the run.

Two further notes from the port:

* **A guard that is false at run time appears to be enough, but is not.** An unassociated
  pointer component referenced under `if (associated(...))` seems to need no map, and it does
  survive in practice. The next section shows why relying on that is a mistake: the standard
  leaves the unmapped case *undefined*, and the guard is in any case a property of the
  *configuration*, not of the code. The `visc%Kv_shear` fault was invisible in one MOM6 test
  case (`ADIABATIC`, so the pointer is never allocated) and fired on every column in another.
  When a new configuration crashes in a region that has been passing for weeks, diff the two
  configurations' parameters and look for a guard that flipped -- bisecting commits will not
  find it, because the bug is present in every commit since the region was offloaded.
* **Mapping the components without the parent worked in every variant tried here**, i.e.
  `map(to: st%alloc, st%ptr)` with no bare `st`. The MOM6 port nevertheless needed the parent
  named as well, so naming both is the habit worth keeping.

### `associated()` on the device does not mean what you think

The bullet above says an unassociated pointer component under `if (associated(...))` needs no
map. That is true, but it is only half the story, and the other half is nastier. The guard is
evaluated *on the device*, against whatever descriptor the device copy of the struct holds —
so the question is whether the device and the host agree. `./null_component` measures it.

Columns are `associated(p)`, `allocated(a)`, `associated(live_p)`, `allocated(live_a)`, where
`p`/`a` are deliberately null/unallocated and `live_p`/`live_a` really do have storage:

| case | mapping | p | a | live_p | live_a |
| :--- | :--- | :-- | :-- | :-- | :-- |
| — | host truth | 0 | 0 | **1** | 1 |
| 1 | `map(to: st)` on the construct | 0 | 0 | **0** | 1 |
| 2 | `enter data map(to: st)` then `map(to: st%p, st%a)` | 0 | 0 | **1** | 1 |
| 3 | `map(to: st, st%p, st%a)` on the construct | 0 | 0 | **0** | 1 |
| 8 | as 1, teams-distribute instead of plain `target` | 0 | 0 | **0** | 1 |

Three things fall out of that table.

**Null and unallocated components are harmless.** Every configuration reports them false, and
naming them in a map clause is not an error. So a guarded map like MOM6's

```fortran
if (associated(visc%Kv_shear)) then
  !$omp target enter data map(to: visc%Kv_shear)
endif
```

is not required for the *map* to be safe. Keep the guard for hygiene, not for survival.

**An associated component that is not named is reported as unassociated.** That is the
`live_p` column, and it is a silent wrong answer, not a fault: case 7 runs a kernel whose
guarded branch should sum to 1296 and gets 0, with no diagnostic of any kind. Nothing in the
run says a branch was skipped.

**Allocatables and pointers differ.** `allocated(live_a)` is true on the device in every case
while `associated(live_p)` is not, so an `allocated()` guard cannot be used as evidence that an
`associated()` guard behaves the same way.

Then the part that makes this genuinely dangerous. An unmapped pointer component can appear to
work perfectly:

```
$ ./null_component 5     # associated pointer, parent mapped, component never named
    returned without faulting, sum =  1.296000000000000E+03      # the correct answer
$ LIBOMPTARGET_INFO=32 ./null_component 5 | grep -c 'Size=512'
0                                                               # never transferred
```

The right answer, and the array was never copied to the device. Case 9 shows what is actually
happening — map the parent, overwrite the host array without any `target update`, then run the
kernel:

```
$ ./null_component 9
    sum                        =  1.296000000000000E+05          # 129600, not 1296
```

The kernel is reading **host memory directly**. On an integrated GPU sharing physical RAM that
often just works, which is exactly why this class of bug survives testing: the same
construction is what faulted as `visc%Kv_shear` in MOM6, and whether you get the right answer,
a wrong answer, or a memory access fault depends on the allocation rather than on the code.

Practical consequence: never treat a passing test as evidence that a component is correctly
mapped. Check `LIBOMPTARGET_INFO=32` for the transfer, or poison the host copy the way case 9
does.

#### What the standard says, and why you should map it anyway

The measurements above are one permitted outcome, not a guarantee. OpenMP 5.2 §5.8.3, on a
component of a mapped derived type that is *not* itself named in a map clause on the construct:

> If it has the POINTER attribute, the map clause treats its association status as if it is
> **undefined**

Undefined -- not disassociated. So the device-side `associated()` that a guarded kernel branch
depends on has no defined value at all when the component is unmapped, and the "device reports
false" seen above is simply what this implementation happens to do. OpenMP 6.0 §7.9.6 replaces
that sentence with "it is **attach-ineligible**", which leaves the device holding a byte copy of
the host descriptor -- a stale host address, which is exactly the case-9 behaviour.

Naming the component is what buys a guarantee. OpenMP 5.2 §13.8 constrains the mapped case:

> If the association status of a list item with the POINTER attribute that appears in a map
> clause on the construct is disassociated upon entry to the target region, the list item must
> be disassociated upon exit from the region.

That is the property a guarded kernel actually wants, and it is only available once the
component appears in a clause. Hence the rule worth following:

**Name every component the region touches, unconditionally. Do not wrap the map in
`if (allocated(...))` or `if (associated(...))`.** Mapping one that has no storage transfers
nothing and is explicitly legislated for; omitting one that does have storage is undefined.
MOM6's `MOM_vert_friction.F90` was changed to do this, deleting eight such guards.

Reads through a missing map can go to host memory and quietly produce the right answer, as
above. **Writes are not so forgiving.** In MOM6, registering the `du_dt_visc` / `du_dt_str`
diagnostics associates `ADp%du_dt_str`, whose components were in no map clause, and the run
dies in the first tridiagonal solve:

```
OFFLOAD ERROR: memory access fault by GPU 1 (agent 0x24190b10) at virtual address
0x2531f000. Reasons: Page not present or supervisor privilege, Write access to a
read-only page
```

Naming the parent and then the six components fixes it. Same missing map, two different
symptoms depending on whether the kernel reads or writes.

### `reduce()` on `do concurrent` is accepted and ignored on the device

`dc_reduce_test.F90`. F2023 lets a `do concurrent` carry a `reduce()` locality specifier, and
upstream MOM6 uses it behind a `DO_LOCALITY(...)` macro. Under
`-fdo-concurrent-to-openmp=device` it compiles, the kernel launches, and the reduction result
never comes back:

| construct | device | expected |
|---|---|---|
| `do concurrent (i=1:n)` writing `s = 42` | **0** | 42 |
| `do concurrent (i=1:n) reduce(+: s)` | **0** | 1024 |
| the same `reduce()` under `-fdo-concurrent-to-openmp=host` | 1024 | 1024 |
| `!$omp target teams distribute parallel do reduction(+: s)` | 1024 | 1024 |

So it is not that `reduce()` is unimplemented -- the host lowering honours it. The device
lowering drops it, and drops plain scalar write-back with it. Two neighbouring traps:

- one file containing both `do concurrent ... reduce(+: int)` and `!$omp ... reduction(+: int)`
  fails to compile with `error: redefinition of symbol named 'add_reduction_i32'`;
- some shapes crash the compiler in `DoConcurrentConversion::genReductions`.

The danger is not a wrong count. In MOM6 the reduced value is `trunc_any`, which decides
whether velocity truncation runs at all, so the silent `.false.` would have disabled it.
**Use `!$omp target teams distribute parallel do ... reduction(...)`.** `reduction(.or.: a, b)`
over two logicals works. If the loop also takes a running `min` into an array element, collapse
only the horizontal indices and leave the vertical one sequential inside the kernel: collapsing
it races, and column-by-column happens to be the host's visit order, so answers stay bitwise
identical.

### The per-thread device stack is 1024 bytes by default

Private *automatic* arrays -- bounds known only at run time, so they can be neither
registerised nor given static storage -- live on the device stack, and ROCm gives each
thread 1 KB of it. Overflow it and the kernel returns wrong numbers. There is no fault, no
warning, no diagnostic of any kind, and the numbers differ from run to run.

`vertvisc_coef` hits this squarely: each of its regions carries six `real, dimension(SZK_(GV))`
automatics in `private()`. `./private_stack` reduces that to six automatics and nothing else,
and runs the same kernel twice:

```
  FAIL  GPU run 1 vs CPU — 4096 of 4096 points differ
  FAIL  GPU run 2 vs CPU — 4096 of 4096 points differ
  FAIL  GPU run 1 vs GPU run 2 — 4096 of 4096 points differ
```

That third line is the diagnostic worth remembering. Two bitwise-identical runs of the same
binary disagreeing with *each other* rules out an arithmetic or porting mistake immediately,
and points at memory rather than at the code.

Sweeping `LIBOMPTARGET_STACK_SIZE` at nz=24, three checks per run, repeated three times with
identical results:

| `LIBOMPTARGET_STACK_SIZE` | result |
| :--- | :--- |
| unset (1024) | silent corruption, differs run to run |
| 2048, 4096, 6144 | silent corruption |
| 8192 | passes |
| 12288, 16384, 24576 | memory access fault |
| 32768, 65536, 131072 | passes |

**The response is not monotonic**, which is worth knowing before you go hunting for a
threshold: 8192 passes, 16384 faults outright, 32768 passes again. The next section works out
why.

The array bytes alone do not predict the requirement either. Six automatics at nz=24 is
1152 bytes, but the kernel needs about 7 KB before it stops corrupting, and at nz=4 it is only
192 bytes and still fails at the 1 KB default. It is not register spill -- the kernel's
metadata reports `.sgpr_spill_count: 0` and `.vgpr_spill_count: 0`. Each runtime-bounded array
simply gets its own aligned allocation in the dynamic stack frame, so the frame is several
times the size of the data in it. That matches the MOM6 experience, where the failing
configuration had nz=2.

### Which values of `LIBOMPTARGET_STACK_SIZE` are actually valid

Three of the four boundaries are exactly derivable; the fourth is not, and that is the useful
thing to know. ROCm ships the plugin source, so these can be read rather than guessed:
`/opt/rocm-7.2.3/lib/llvm/share/gdb/python/ompd/src/offload/plugins-nextgen/`.

**It does nothing at all for most kernels.** `rtl.cpp` sets the dispatch packet's private
segment to `Kernel.usesDynamicStack() ? max(getPrivateSize(), StackSize) : getPrivateSize()`.
Only kernels whose metadata says `.uses_dynamic_stack: true` are affected -- runtime-bounded
private arrays are what sets that flag. Every one of the mapping and `do concurrent`
reproducers here is immune to the variable no matter what it is set to.

**The default is 1024 bytes**, hardcoded as `uint32_t StackSize = 1024 /* 1 KB */`
(`rtl.cpp:4991`), with the comment "in conformity to hipLimitStackSize". Note that upstream
LLVM uses `16 * 1024` instead -- so a value that is fine on a stock LLVM build may corrupt on
ROCm, and the same source read on the wrong branch will mislead you.

**Parsing is `istringstream >> uint64_t`**: plain decimal bytes, no units and no hex, so
`0x8000` parses as 0 and an unparseable value is silently ignored, leaving the 1 KB default.

**Granularity is 8 bytes.** ROCr rounds the per-work-item size up to `(64*4)/WavefrontSize`,
which is 8 B on gfx11 wave32 and 16 B on CDNA wave64.

**The upper bound is 262112, not the 262136 the plugin advertises.** `MaxThreadScratchSize` is
`((64*4)/WavefrontSize) * ((1<<15)-1)` = 262136 for gfx11 wave32 (`((256*4)/WavefrontSize) *
((1<<13)-1)` = 131056 on gfx90a/gfx942). Anything larger is clamped to it with an "AMDGPU
message: Scratch memory size will be set to ..." on stderr. But ROCr's own `MAX_WAVE_SCRATCH`
works out to 262112 per work-item, 24 bytes lower, and it refuses rather than clamps:

| value | result |
| :--- | :--- |
| 262112 | passes -- the true ceiling |
| 262120, 262136 | `HSA_STATUS_ERROR_OUT_OF_RESOURCES`, fatal |
| 262144 | clamped to 262136, then fatal anyway |

So the plugin's own cap is unusable: the band between the two limits is accepted and then dies.

**The fault band is the scratch-reclaim threshold.** ROCr commits scratch for the whole device
at full occupancy -- here 12 CUs x 32 waves/CU x 32 lanes = **12288 work-items** (from
`rocminfo`) -- and if that commitment exceeds `HSA_SCRATCH_SINGLE_LIMIT` it switches from
queue-retained scratch to the "large"/single-use path, which faults on this gfx1103. So:

```
   AlignUp(12288 * RoundUp8(value), 65536)  >  146800640   ->  fault
```

That predicts a largest safe value of exactly **11944**, and it was confirmed out of sample:

| value | predicted | observed |
| ---: | :--- | :--- |
| 11936, 11940, 11944 | pass | pass |
| 11945, 11946, 11952 | fault | fault |

The mechanism is causally established rather than merely fitted -- moving
`HSA_SCRATCH_SINGLE_LIMIT` moves the band with it. Dropping it to 8 MB drags the fault down
over 8192, and raising it to 1 GB (or setting `HSA_NO_SCRATCH_RECLAIM=1`, which forces
`large = false`) clears 16384/24576/30720.

**What is not derivable:** why values above about 32000 work again, and why isolated values
inside the safe window still corrupt. At nz=24, 7936, 8320, 8576, 9216, 9472, 9728 and 10624
all corrupt reproducibly while their immediate neighbours 128 bytes either side are fine, and
all of them pass at nz=4. Neither scratch knob fixes them. They are demand-dependent, so they
behave like the runtime delivering less usable stack than asked for, by an amount that is not
monotonic in the request.

**The practical rule**, then:

* Stay at or below **11944** to keep ROCr on the retained-scratch path, and above whatever your
  kernel actually needs -- for the six-array reproducer that is about 7 KB at nz=24.
* Do not treat a bigger number as a safer number. Above 11944 you are relying on the
  single-use path, which faults for most of 12032-31744 and only works again by accident.
* Verify the specific value with your own kernel at your own `nz`, because of the isolated bad
  values. `./private_stack <nz>` is the check.
* **Do not sit exactly on the boundary.** 11944 commits precisely the 140 MB limit, and MOM6's
  `benchmark` never got past its first timestep in 25 minutes there, while the same run takes
  134 s at 65536.

#### What ROCm 10.0.0 changes

Rechecked against a ROCm 10.0.0 userspace install (`amdflang 23.0.0git`,
`libomptarget.so.23.0git`) on the same gfx1103, reproducer rebuilt with its compiler.

**The 1024-byte default is unchanged.** ROCm 10 ships no plugin source, but the constructor
still stores it: `movl $0x400, 0xb90(%rbx)` inside
`AMDGPUDeviceTy::AMDGPUDeviceTy(...)`, and `getDeviceStackSize` reads that same offset with a
32-bit load, so it is still a `uint32_t` initialised to 1 KB. Confirmed at runtime too --
unset and an explicit 1024 corrupt identically, and 1280 is the first value that passes.

**The fault band is gone**, and genuinely fixed rather than moved: every value from 11952
through 32768 now passes, and forcing `HSA_SCRATCH_SINGLE_LIMIT=8388608` -- which on 7.2.3
dragged the fault down over 8192 -- leaves everything passing. The single-use scratch path
appears to be repaired.

**The ceiling bug survives.** 262112 passes; 262120 and 262136 still die with
`HSA_STATUS_ERROR_OUT_OF_RESOURCES`; anything larger is still clamped to 262136 and then dies
anyway. The 24-byte disagreement between the plugin's cap and ROCr's real limit is unchanged.

**The requirement dropped sharply**, which is a codegen change rather than a runtime one: the
same source at nz=24 needs about 7 KB built with ROCm 7.2.3 and only 1280 bytes built with
ROCm 10. The dynamic stack frame for runtime-bounded arrays got much tighter. Note that 1024
still fails even so -- the default remains too small for a kernel with six small automatics.

So on ROCm 10 the practical rule collapses to the two derivable bounds: above what the kernel
needs, and at or below 262112. The awkward middle disappears.

#### Staying on 7.2.3

For anything with a realistic column count the window is simply too small. MOM6's
`vertvisc_coef` at nz=22 produces wrong answers (energy 1.98 against 0.547, 313 velocity
truncations) at 11264, so it cannot be run on the retained path at all with the stock limit.
The principled fix is to move the band rather than step over it:

```sh
export HSA_SCRATCH_SINGLE_LIMIT=1073741824   # 1 GB, up from the 140 MB default
export LIBOMPTARGET_STACK_SIZE=65536
```

which puts a 65536-byte stack back on the retained path (12288 x 65536 = 805 MB < 1 GB) and
reproduces MOM6 `benchmark` bitwise. The cost is that the scratch is committed device-wide:
12288x the per-work-item value on this GPU, so 65536 commits 805 MB of the shared system RAM
an APU uses as VRAM.

The general diagnostic for all three of these: re-run the *same binary* with
`OMP_TARGET_OFFLOAD=DISABLED` and `OMP_NUM_THREADS=1`. A bitwise pass there proves the
arithmetic is right and the fault is offload-specific. Host-*threaded* is not a valid oracle
-- it has its own failures -- so use serial.

## Summary

* `amdflang` supports basic OpenMP offload, and to get good results, the long form 
  directives must be used and `loop` and `do concurrent` must be avoided. This is very 
  unfortunate, as `nvfortran` prefers `loop` and `do concurrent` over the long form OpenMP 
  directive.
* `amdflang` doesn't handle nested parallelism as well as `nvfortran` does.
  However, it seems that the most recent ROCm (7.2.3) can get ok performance.
* `amdflang` benefits more strongly from explicitly setting number of teams and threads.
  This doesn't present as an obstacle, as it seems that it also benefits `nvfortran`.
* The asynchronous default behaviour of OpenMP kernel launches in `amdflang` means
  `!$omp taskwait` must be more diligently used, or the `OMPX_FORCE_SYNC_REGIONS`
  environment variable must be set.
* Despite AMD OpenMP kernels being async, we can't expect them to be able to hide the
  kernel launch overhead for many small kernels.
* All of the correctness traps above fail *silently*. Three of them (the `exit`
  lowering, the device stack, and an unmapped pointer component read out of host
  memory) can produce a bitwise match against the CPU baseline while being completely
  wrong, so a passing regression test is not by itself evidence that a loop was ported.
* `associated()` evaluated inside a kernel is not a reliable mirror of the host. An
  associated pointer component that is not named in a map clause is reported as
  *unassociated* on the device, so a guarded branch is skipped and the kernel quietly
  computes something else. `allocated()` on an allocatable component does not behave
  this way, so one cannot be used to reason about the other.
* A missing map shows up as a silent wrong answer when the kernel *reads* the component
  and as a hard `memory access fault ... Write access to a read-only page` when it
  *writes* one. Only the second is self-announcing, and only if the run reaches it.
* No scalar written inside a `do concurrent` survives the device lowering, `reduce()`
  included -- it reads back as its pre-loop value. That makes `do concurrent` unusable
  for reductions here, which is a further reason the long-form directives are the ones
  to write on AMD.

The main implications for MOM6 porting are:
* The latest ROCm available should be used. At the time of testing, Pawsey's newest ROCm
  module was 6.4.1.
* Since NVIDIA does better with `do concurrent` and OpenMP `loop`, but AMD
  does better with `distribute parallel do`, there must be some conversion between
  them (maybe macros).
* Every ported region needs a positive check that a kernel was actually dispatched
  (`LIBOMPTARGET_INFO=16`), not just a checksum comparison.
* Nothing in the NVIDIA port's map clauses can be trusted to carry over, because
  `mem:managed` means it largely doesn't have any. Every derived-type component a region
  touches -- including the ones reached several call levels down, and including the ones
  behind guards that happen to be false in the configuration being tested -- has to be
  mapped explicitly.
* `LIBOMPTARGET_STACK_SIZE` needs to be set for every run, and belongs in whatever wrapper
  script is used to launch the model rather than being left to the user to remember. Set it to
  a validated value at or below 11944 on this device, not to the largest number that happens
  to work -- above that threshold ROCr switches scratch strategy and most of the range faults.
  A deep ocean configuration may need more than 11944, in which case raise
  `HSA_SCRATCH_SINGLE_LIMIT` rather than stepping over the band.

## Bandwidth microbenchmark

`bandwidth_test.F90` measures the raw cost of `map(to:)` / `map(from:)` / `map(tofrom:)`
data transfers on the Radeon 780M iGPU (`gfx1103`), separately from any compute. It sweeps
array sizes from 4 KiB to 256 MiB through five patterns -- cold `map(to)`, cold `map(from)`,
cold `map(tofrom)`, and warm `target update to`/`from` against a device buffer that's
allocated once outside the timed loop -- plus a kernel-launch-overhead floor (a single
already-resident element, no bulk transfer). Bandwidths are best-of-5, in GB/s (10^9 B/s);
`map(tofrom)`'s GB/s is computed against `2*bytes` since it crosses the bus both ways. Build
and run it the same way as the other tests here:

```sh
make bandwidth FC=/opt/rocm-7.2.3/bin/amdflang \
  FFLAGS="-fopenmp --offload-arch=gfx1103 -O3 -ffp-contract=off -fdo-concurrent-to-openmp=device"
./bandwidth
```

(`amdflang` isn't on `PATH` by default on this machine -- pass the full path as above, or add
`/opt/rocm-7.2.3/bin` to `PATH` first.)

Measured on the Radeon 780M, bare metal, ROCm 7.2.3:

```
Kernel launch overhead (single resident element, best of 5 runs):     12.159 us

       bytes    map(to)   map(from)  map(tofrom)  update(to) update(from)
        4096       0.158       0.158       0.222       0.051       0.100
       16384       0.608       0.608       0.799       0.205       0.390
       65536       1.883       2.051       2.581       0.818       1.462
      262144       5.848       6.247       7.505       2.916       4.780
     1048576      12.785      16.916      17.047       8.267      13.450
     4194304       7.269       8.290      11.748      16.644      23.425
    16777216      10.434      10.552      14.457      21.619      23.140
    67108864      10.994      10.414      13.729      19.785      18.143
   268435456      10.852      10.404      13.464      23.974      22.978
```

Small transfers are dominated by the ~12 us launch floor, not the interconnect -- a 4 KiB
transfer at "0.16 GB/s" is really just ~12 us of fixed overhead moving 4 KiB, not a slow
copy. Bandwidth climbs through the tens of KiB-MiB range and settles around 10-14 GB/s
for the cold `map(...)` patterns and 18-24 GB/s for the warm `target update` patterns once
allocation overhead is out of the timed region -- both well under the theoretical DDR
bandwidth of this APU's shared memory, since these transfers go through the OpenMP runtime's
own copy path rather than a raw `memcpy`. The dip at 4 MiB across every pattern (e.g.
`map(to)` falling from 12.8 to 7.3 GB/s) repeats across runs and is most likely a cache/TLB
transition rather than noise, but wasn't investigated further.

### Two more traps, found writing this benchmark

Neither is about wrong answers -- one is a hard crash, the other is a silently wrong
diagnostic -- but both are the same shape as the traps above: nothing in the build log says
anything is wrong, and the code compiles clean.

**A `!$omp target` (or `target teams distribute parallel do`) construct written directly
inside a host `do` loop can compile fine and then segfault on every call**, with:

```
"PluginInterface" error: Failure to init kernel: "requested object was not found in the
binary image" ... HSA_STATUS_ERROR_INVALID_SYMBOL_NAME: There is no symbol with the given name.
omptarget error: Failed to load kernel __omp_offloading_...
```

The host-side dispatch call survives compilation, but the matching kernel symbol is missing
from the device code object -- something in the device linking path drops it when the
construct's enclosing scope is a runtime-trip-count host loop rather than a plain subroutine
body. The fix is structural, not a flag: put the `target` construct in its own subroutine and
call that subroutine from the timing loop one level up, exactly like every other test in this
directory already does (`run_av_rem_omp` etc. in `av_rem_test.F90` are called from a
`do irun = 1, n_runs` loop in the *caller*; the directive itself is never lexically inside
that loop). First-timed variants of every one of the five patterns above hit this before being
restructured this way; none did afterward.

**`omp_lib`'s `omp_is_initial_device()` reports HOST from inside a `target teams distribute
parallel do`, even while `LIBOMPTARGET_INFO=16` shows the kernel actually launching on the
device.** Rerunning with that env var set showed real `Launching kernel ... AMDGPU device 0`
lines for the exact call that had just reported itself as running on the host. `dc_exit_test.F90`
in this same directory already carries the workaround, for a different reason (getting past
`do concurrent`'s purity rule) -- it declares its own C-binding interface straight to the
`omp_is_initial_device` symbol instead of using the one `omp_lib` provides:

```fortran
interface
  integer(c_int) function is_initial_device() bind(C, name="omp_is_initial_device")
    import :: c_int
  end function is_initial_device
end interface
```

That version reports correctly here too. So: don't trust `omp_lib`'s device-dispatch check
inside a `teams`/`parallel` nesting on this compiler -- use the raw C binding, and treat a
passing dispatch check done the `omp_lib` way as inconclusive rather than as proof of a host
fallback.
