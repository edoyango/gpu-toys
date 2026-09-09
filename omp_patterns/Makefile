FC     ?= amdflang
FFLAGS ?= -fopenmp --offload-arch=native -O3 -ffp-contract=off  -fdo-concurrent-to-openmp=device

UTILS = test_utils.F90

.PHONY: all clean

all: av_rem col_norm repro_tests dc_exit map_components private_stack null_component \
     dc_reduce bandwidth

av_rem: $(UTILS) av_rem_test.F90
	$(FC) $(FFLAGS) -o $@ $^

col_norm: $(UTILS) col_norm_test.F90
	$(FC) $(FFLAGS) -o $@ $^

repro_tests: $(UTILS) repro_tests.F90
	$(FC) $(FFLAGS) -o $@ $^

dc_exit: $(UTILS) dc_exit_test.F90
	$(FC) $(FFLAGS) -o $@ $^

map_components: $(UTILS) map_components_test.F90
	$(FC) $(FFLAGS) -o $@ $^

private_stack: $(UTILS) private_stack_test.F90
	$(FC) $(FFLAGS) -o $@ $^

null_component: $(UTILS) null_component_test.F90
	$(FC) $(FFLAGS) -o $@ $^

dc_reduce: dc_reduce_test.F90
	$(FC) $(FFLAGS) -o $@ $^

bandwidth: $(UTILS) bandwidth_test.F90
	$(FC) $(FFLAGS) -o $@ $^

clean:
	rm -f av_rem col_norm repro_tests dc_exit map_components private_stack \
	      null_component dc_reduce bandwidth *.mod
