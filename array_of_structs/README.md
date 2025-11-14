MOM6 tracer advection uses array of structs with array members and I was noticing pretty poor performance.
This demonstrates that `nvfortran` struggles to compile optimal code when using such a data structure.

## TL;DR

Don't reference arrays of structs with dynamic array members inside OpenMP offloaded loops.

## nsys output

Each of the loops are wrapped in `nvtx` markers to highlight the differences. The times shown include
data transfers upon entry/exit into each `omp target`/`do concurrent` loop.
```
 Time (%)  Total Time (ns)  Instances     Avg (ns)         Med (ns)        Min (ns)       Max (ns)     StdDev (ns)   Style       Range     
 --------  ---------------  ---------  ---------------  ---------------  -------------  -------------  -----------  -------  --------------
     59.3    1,914,988,980          1  1,914,988,980.0  1,914,988,980.0  1,914,988,980  1,914,988,980          0.0  PushPop  :kjmi         
     34.3    1,107,055,904          1  1,107,055,904.0  1,107,055,904.0  1,107,055,904  1,107,055,904          0.0  PushPop  :mkji flat    
      6.4      206,199,920          1    206,199,920.0    206,199,920.0    206,199,920    206,199,920          0.0  PushPop  :mkji with ptr
```
Considering purely kernel compute time: 
```
 Time (%)  Total Time (ns)  Instances     Avg (ns)         Med (ns)        Min (ns)       Max (ns)     StdDev (ns)           Name          
 --------  ---------------  ---------  ---------------  ---------------  -------------  -------------  -----------  -----------------------
     58.1    1,539,002,761         40     38,475,069.0     37,879,328.5     29,895,548     50,929,747  3,241,844.6  nvkernel_MAIN__F1L37_2_
     41.8    1,106,699,922          1  1,106,699,922.0  1,106,699,922.0  1,106,699,922  1,106,699,922          0.0  test_51_gpu            
      0.1        2,041,848        400          5,104.6          4,384.0          4,319        103,039      7,326.5  nvkernel_MAIN__F1L64_4_
```
The worst case is the `kjmi` loop, which is a rough sketch of the loop ordering seen in MOM6 tracer advection.
The best case is the `mkji with ptr` loop. The fact that `mkji flat` loop is slow probably indicates that the
use of arrays of structs with allocatable array members is a terrible idea.

<img width="926" height="159" alt="image" src="https://github.com/user-attachments/assets/6a302baa-64c1-4348-8b64-61dd75af73b6" />

When looking at the nsys profiling data, you can see page fault data transfers happening within the kernel.
