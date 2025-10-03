FC     ?= gfortran
FFLAGS ?= -fopenmp -O2

SRC = src/cp2kmd2g96.f90

TARGET = cp2kmd2g96.x

all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -rf $(TARGET)

.PHONY: all clean
