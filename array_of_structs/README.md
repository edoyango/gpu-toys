MOM6 tracer advection uses array of structs with array members and I was noticing pretty poor performance.

~~This demonstrates that `nvfortran` struggles to compile optimal code when using such a data structure.~~

## TL;DR

~~Don't reference arrays of structs with dynamic array members inside OpenMP offloaded loops.~~

It turns out you have to make sure you transfer both the array of structs and the arras that are members of the structs.

## nsys output

Each of the loops are wrapped in `nvtx` markers to highlight the differences. The times exclude any data transfer time
at setup/teardown.
```
 Time (%)  Total Time (ns)  Instances    Avg (ns)      Med (ns)     Min (ns)    Max (ns)   StdDev (ns)   Style       Range     
 --------  ---------------  ---------  ------------  ------------  ----------  ----------  -----------  -------  --------------
     38.5       12,264,220          1  12,264,220.0  12,264,220.0  12,264,220  12,264,220          0.0  PushPop  :kjmi         
     33.8       10,789,946          1  10,789,946.0  10,789,946.0  10,789,946  10,789,946          0.0  PushPop  :mkji with ptr
     27.7        8,840,687          1   8,840,687.0   8,840,687.0   8,840,687   8,840,687          0.0  PushPop  :mkji flat
```
The best case is the `kjmi flat` loop, and the worst case is the `kmji` loop - the latter being a rough sketch of what
happens in MOM6 tracer advection. Having a flatter loop is clearly beneficial.

## without proper data mapping

Considering purely CUDA kernel time: 
```
 Time (%)  Total Time (ns)  Instances     Avg (ns)         Med (ns)        Min (ns)       Max (ns)     StdDev (ns)           Name          
 --------  ---------------  ---------  ---------------  ---------------  -------------  -------------  -----------  -----------------------
     58.1    1,539,002,761         40     38,475,069.0     37,879,328.5     29,895,548     50,929,747  3,241,844.6  nvkernel_MAIN__F1L37_2_
     41.8    1,106,699,922          1  1,106,699,922.0  1,106,699,922.0  1,106,699,922  1,106,699,922          0.0  test_51_gpu            
      0.1        2,041,848        400          5,104.6          4,384.0          4,319        103,039      7,326.5  nvkernel_MAIN__F1L64_4_
```
The worst case is still the `kjmi` loop, but now takes over 100x longer. Looking at the nsys profiling output reveals that
data transfers are happening due to page faults during kernel runtime causing data transfers to occur mid-kernel.

<img width="926" height="159" alt="image" src="https://github.com/user-attachments/assets/6a302baa-64c1-4348-8b64-61dd75af73b6" />
