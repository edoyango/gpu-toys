# OpenMP loop patterns on AMD and NVIDIA

This subdirectory contains tests for OpenMP + do concurrent loop patterns for both AMD and NVIDIA.
The reason for this comparison is that we've ported MOM6 code using OpenMP + do concurrent for
NVIDIA only, and I would like to get a sense of how well our code would translate to AMD. Hence,
the tests are heavily reduced versions of some common loop patterns used so far.

Hardware:

* AMD: MI250x single GCD using ROCm 6.4.1.
  * Latest ROCm is 7.2.3, but the system used doesn't have newer ROCm.
  * Code compiled with `amdflang -fopenmp --offload-arch=gfx90a -Mnofma -O3`
  * 28 TFLOPS, 1.6TB/s
* NVIDIA: V100 SXM using NVHPC 25.9.
  * Code compiled with `nvfortran -mp=gpu -gpu=mem:separate -Mnofma -Mnovect -O4 -Minline=flux_elem`
  * 7.5 TFLOPS, 0.9TB/s

The `-ffp-contract=off` and `-Mnofma -Mnovect` flags are to ensure that results are bitwise
identical with CPU for NVIDIA. Interestingly, excluding these flags *slows down* the NVIDIA GPU code.
`-Minline=flux_elem` is included in the NVIDIA compilation as `amdflang` automatically inlines at `-O3`
(can be checked with `-Rpass=inline`). Furthermore, there is a compiler bug affecting some of these
tests which require inlining functions to get correct results.

It's important to note that MI250x is a few years newer than the V100, and a more appropriate
comparison would be with the A100s. But the latter were harder to get a hold of at the time.

All timings show best of 5 runs using 512x512x100 grid.

Note that wherever the `loop` construct is used in the below examples, the AMD code actually uses
the older `teams distribute parallel do`. On the otherhand, NVIDIA uses `loop`. There are two
reasons for this difference:
1. The version of ROCm used here gives the wrong answers when using `loop`
  * When tested on my laptop with a newer ROCm, `loop` does work correctly.
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
| MI250x ROCm 6.4.1 | 2.208     | 2.453     | 2.507       | 4.391       | 15.831      |
| V100 NVHPC 25.9   | 1.168     | 1.454     | 1.658       | 4.089       | 12.060      |

Despite being older and weaker, the V100 beats the MI250x in all problem sizes. 

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
| MI250x ROCm 6.4.1 | 0.548     | 0.697     | 0.722       | 1.496       | 6.671       |
| V100 NVHPC 25.9   | 0.374     | 0.667     | 0.859       | 2.234       | 8.881       |

And in both cases, the compiler does a better job of optimising. Furthermore, the AMD setup is now
doing better than the tested NVIDIA setup by a meaningful margin at larger problem sizes.

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

|                   | 32x32x100 | 64x64x100 | 128x128x100 | 256x256x100 | 512x512x100 |
| :---              | ---:      | ---:      | ---:        | ---:        | ---:        |
| MI250x ROCm 6.4.1 | 0.956     | 1.823     | 3.713       | 7.440       | 16.536      |
| V100 NVHPC 25.9   | 0.038     | 0.058     | 0.071       | 0.169       | 0.681       |

The AMD setup is clearly struggling here, being at least an order of magnitude slower than the
NVIDIA setup.

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
| V100 NVHPC 25.9   | 0.014     | 0.016     | 0.041       | 0.140       | 0.541       |

The timings look like the difference between loop pattern 1a and 1b. 
