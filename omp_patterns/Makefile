FC     ?= amdflang
FFLAGS ?= -fopenmp -offload-arch=native -O3 -ffp-contract=off

UTILS = test_utils.F90

.PHONY: all clean

all: av_rem col_norm repro_tests

av_rem: $(UTILS) av_rem_test.F90
	$(FC) $(FFLAGS) -o $@ $^

col_norm: $(UTILS) col_norm_test.F90
	$(FC) $(FFLAGS) -o $@ $^

repro_tests: $(UTILS) repro_tests.F90
	$(FC) $(FFLAGS) -o $@ $^

clean:
	rm -f av_rem col_norm repro_tests *.mod
