# OpenMP loop patterns on AMD and NVIDIA

This subdirectory contains tests for OpenMP + do concurrent loop patterns for both AMD and NVIDIA.
The reason for this comparison is that we've ported MOM6 code using OpenMP + do concurrent for
NVIDIA only, and I would like to get a sense of how well our code would translate to AMD. Hence,
the tests are heavily reduced versions of some common loop patterns used so far.

Setups:

* AMD: MI250x single GCD using ROCm 7.2.3.
  * This is run with singularity - ROCm 6.4.1 and kernel version 6.4 on host.
  * Code compiled with `amdflang -fopenmp --offload-arch=gfx90a -Mnofma -O3 -ffp-contract=off -fdo-concurrent-to-openmp=device`
  * 28 TFLOPS, 1.6TB/s
* NVIDIA: V100 SXM using NVHPC 25.9.
  * Code compiled with `nvfortran -mp=gpu -gpu=mem:separate -Mnofma -Mnovect -O4 -Minline=flux_elem -stdpar=gpu`
  * 7.5 TFLOPS, 0.9TB/s

The `-ffp-contract=off` and `-Mnofma -Mnovect` flags are to ensure that results are bitwise
identical with CPU. Interestingly, excluding these flags *slows down* the NVIDIA GPU code.
`-Minline=flux_elem` is included in the NVIDIA compilation as `amdflang` automatically inlines at `-O3`
(can be checked with `-Rpass=inline`). Furthermore, there is a compiler bug affecting some of these
tests which require inlining functions to get correct results.

It's important to note that MI250x is a few years newer than the V100, and a more appropriate
comparison would be with the A100s. But the latter were harder to get a hold of at the time.

All timings show best of 5 runs using 512x512x100 grid.

Note that wherever the `loop` construct is used in the below examples, the AMD code actually uses
the older `teams distribute parallel do`. On the otherhand, NVIDIA uses `loop`. There are two
reasons for this difference:
1. with `amdflang`, `loop` didn't always give the right answers or they were slow.
2. NVIDIA performs significantly faster with `loop` than the older construct.

## Loop pattern 1a: Outer target teams with parallel inner ij loops + serial loops

This is probably the most complex loop pattern used so far. The `target teams` region is started
and then an initialisation ij loop is run. Following is a newton-raphson iterative loop that ought
to be serialised for every thread. Inside the newton-raphson loop is more parallel ij loops, with
one of them being inside a serialised k loop. The k loop needs to be serialised because there are
sum reductions and bitwise reprodicubility is a goal.

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

Despite being older and weaker, the V100 beats the MI250x in all problem sizes. 

In this scenario, I couldn't get the `loop` construct working on the inner loops in the AMD setup.

Worth noting that when using ROCm 6.4.1, the AMD setup was about 10-15% slower than the above reported
times.

## Loop pattern 1b: Outer ij loop with inner serial loops

The function of this version of the code is the same as above, except structured in a more GPU-friendly
way. Instead of the ij loops being insde, the ij loop is outside with inner loops only being the serial
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

This loop pattern is related to, but simpler than the first loop pattern. This simple reduces a 3d
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

Notably, **when using ROCm 6.4.1, these timings were an order magnitude worse then 7.2.3**.
Which highlights how the AMD compilers are improving. 

Even with the newer ROCm, the AMD setup is slower than NVIDIA by approx 2x - except for the
256x256x100 test. This might be because the block size used in AMD is 256, which matches up 
nicely with that problem size.

I did notice that the 512x512x100 was reduced to ~0.8ms if number of threads was set to match
the `nx` dimension. In contrast, the 256x256x100 time doubled, and the smaller sizes
improved a negligable amount. Now, given that problem sizes are unlikely to line up with
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
support ROCm, but ROCm 7+ does.

For the jki loop:

|                   | 32x32x100 | 64x64x100 | 128x128x100 | 256x256x100 | 512x512x100 |
| :---              | ---:      | ---:      | ---:        | ---:        | ---:        |
| MI250x ROCm 7.2.3 | 1.433     | 2.770     | 5.247       | 10.704      | 23.818      |
| V100 NVHPC 25.9   | 0.019     | 0.021     | 0.045       | 0.140       | 0.525       |

Which is clearly awful on the AMD setup. By comparison, the NVIDIA setup is close to
the OpenMP times. When inspecting the kernel launches on the AMD setup with the
environment varialbe `LIBOMPTARGET_KERNEL_TRACE=1`, the output shows that the kernel is
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
  `!$omp taskwait` much be more diligently used, or the `OMPX_FORCE_SYNC_REGIONS`
  environment variable must set.
* Despite AMD OpenMP kernels being async, we can't expect them to be able to hide the
  kernel launch overhead for many small kernels.

The main implications for MOM6 porting are:
* latest ROCm should be used (Pawsey has 6.4.1 as the latest ROCm as a module)
* Since NVIDIA does better with `do concurrent` and OpenMP `loop`, but AMD
  does better with `distribute parallel do`, there must be some conversion between
  them (maybe macros).

