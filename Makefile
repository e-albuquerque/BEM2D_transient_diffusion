# Makefile for Modernized 2D Transient Diffusion BEM (DCC)

FC = gfortran
FFLAGS = -fopenmp -O3 -fdefault-real-8 -Wall -Wextra
LIBS = -llapack -lblas

ifeq ($(OS),Windows_NT)
    TARGET = dcc_modern.exe
    LDFLAGS = -fopenmp -static
    RM = cmd /C del /Q /F
    PATH := C:\msys64\ucrt64\bin;$(PATH)
else
    TARGET = dcc_modern
    LDFLAGS = -fopenmp
    RM = rm -f
endif

OBJS = data_mod.o solver_mod.o assembly_mod.o io_mod.o main.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

data_mod.o: data_mod.f90
	$(FC) $(FFLAGS) -c $< -o $@

solver_mod.o: solver_mod.f90 data_mod.o
	$(FC) $(FFLAGS) -c $< -o $@

assembly_mod.o: assembly_mod.f90 data_mod.o solver_mod.o
	$(FC) $(FFLAGS) -c $< -o $@

io_mod.o: io_mod.f90 data_mod.o
	$(FC) $(FFLAGS) -c $< -o $@

main.o: main.f90 data_mod.o assembly_mod.o io_mod.o
	$(FC) $(FFLAGS) -c $< -o $@

clean:
ifeq ($(OS),Windows_NT)
	-cmd /C del /Q /F *.o *.mod $(TARGET) 2>nul
else
	$(RM) *.o *.mod $(TARGET)
endif

.PHONY: all clean
