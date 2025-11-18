This test replicates the code path followed in MOM6 horizontal viscosity for the double gyre case.
In its current form, it performs relatively poorly compared to other components when ported to the
GPU. 

The main caveat is that this is heavily trimmed down from the original horizontal viscosity code. 
Optimisations that could be readily be applied to the code in this directory could be more
difficult to apply in practice due to the other code not included here (Open Boundary Conditions,
for example).

Nonetheless, it's a useful piece of code to explore because of its scanning like behaviour in the
vertical direction - performing updates to 2d arrays as the loop traverses across the vertical
direction. This scanning behaviour is very common in MOM6 and is used in almost every other major
component e.g. Continuity, Vertical Viscosity, and Barotropic (to a lesser degree).

## Run the test

To run the test, you must have `nvfortran` available in your path and an NVIDIA GPU. 

Run the test with:

```bash
make test
```

The script will compile and run the test code with a 600 x 600 x 75 grid. The elapsed time of the
Horizontal Viscosity kernel will be printed, alongside whether the result matches with the
hardcoded checksum.

```
=== Running GPU test (600x600x75) ===
 Grid dimensions: ni=          600 nj=          600 nk=           75
 Running horizontal_viscosity subroutine...
 
 === RESULTS ===
 Elapsed time:    0.9064581000000000       seconds
 Checksum diffu:      -3484028734114474076
 Checksum diffv:      -8381563518403603703
 
 Sample diffu values:
   diffu(1,1,1) =    8.9382268913856523E-010
   diffu(mid)   =   -1.6015423332802365E-011
 Sample diffv values:
   diffv(1,1,1) =    6.8536408560045682E-011
   diffv(mid)   =    2.1779787969607433E-010

=== Verifying checksums ===
Reference: diffu=-3484028734114474076, diffv=-8381563518403603703
Actual:    diffu=-3484028734114474076, diffv=-8381563518403603703
✓ PASS: Checksums match reference!
```