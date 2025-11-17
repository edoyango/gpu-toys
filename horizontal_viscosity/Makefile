# Makefile for horizontal viscosity GPU test with checksum validation

# Compiler
FC = nvfortran

# Source file
SRC = main.F90

# Output executable
EXE = test_horvisc_gpu

# Compiler flags
FFLAGS = -mp=gpu -stdpar=gpu -gpu=mem:separate -O4 -Mnofma -Mnovect

# Grid size
NI = 600
NJ = 600
NK = 75

# Reference checksums (update these with your expected values)
REFERENCE_DIFFU ?= -3484028734114474076
REFERENCE_DIFFV ?= -8381563518403603703

# Default target
all: $(EXE)

# Build GPU version
$(EXE): $(SRC)
	$(FC) $(FFLAGS) -o $(EXE) $(SRC)
	@echo "GPU version compiled successfully: $(EXE)"

# Run test and verify checksums
test: $(EXE)
	@echo "=== Running GPU test ($(NI)x$(NJ)x$(NK)) ==="
	@./$(EXE) $(NI) $(NJ) $(NK) > test_output.txt
	@cat test_output.txt
	@echo ""
	@echo "=== Verifying checksums ==="
	@DIFFU=$$(grep "Checksum diffu:" test_output.txt | awk '{print $$3}'); \
	DIFFV=$$(grep "Checksum diffv:" test_output.txt | awk '{print $$3}'); \
	echo "Reference: diffu=$(REFERENCE_DIFFU), diffv=$(REFERENCE_DIFFV)"; \
	echo "Actual:    diffu=$$DIFFU, diffv=$$DIFFV"; \
	if [ "$(REFERENCE_DIFFU)" = "0" ] || [ "$(REFERENCE_DIFFV)" = "0" ]; then \
		echo ""; \
		echo "WARNING: Reference checksums not set!"; \
		echo "Update the Makefile with:"; \
		echo "  REFERENCE_DIFFU ?= $$DIFFU"; \
		echo "  REFERENCE_DIFFV ?= $$DIFFV"; \
		echo "Or override: make test REFERENCE_DIFFU=$$DIFFU REFERENCE_DIFFV=$$DIFFV"; \
		exit 1; \
	fi; \
	if [ "$$DIFFU" = "$(REFERENCE_DIFFU)" ] && [ "$$DIFFV" = "$(REFERENCE_DIFFV)" ]; then \
		echo "✓ PASS: Checksums match reference!"; \
	else \
		echo "✗ FAIL: Checksums do not match reference!"; \
		exit 1; \
	fi

# Clean build artifacts
clean:
	rm -f $(EXE) *.mod *.o test_output.txt

# Help
help:
	@echo "Makefile for horizontal viscosity GPU test"
	@echo ""
	@echo "Targets:"
	@echo "  all   - Build GPU version (default)"
	@echo "  test  - Run test with $(NI)x$(NJ)x$(NK) grid and verify checksums"
	@echo "  clean - Remove compiled executable and test output"
	@echo "  help  - Display this help message"
	@echo ""
	@echo "Reference checksums:"
	@echo "  REFERENCE_DIFFU = $(REFERENCE_DIFFU)"
	@echo "  REFERENCE_DIFFV = $(REFERENCE_DIFFV)"
	@echo ""
	@echo "Override checksums:"
	@echo "  make test REFERENCE_DIFFU=12345 REFERENCE_DIFFV=67890"
	@echo ""
	@echo "Usage:"
	@echo "  make        # Build"
	@echo "  make test   # Build and run test"
	@echo "  make clean  # Clean"

.PHONY: all test clean help