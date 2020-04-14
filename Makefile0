TARGET = main
#LIBS =  -lm  -lgfortran -openmp #-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm -fopenmp


UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
#LIBS = -ldl -lm  -fopenmp 
#CC2 = g++ -fopenmp -g
CC2 = icc -qopenmp -g
#LIBS = -lmkl_intel_lp64  -lmkl_core -ldl -lpthread -lm -qopenmp -llapack -lblas -lmkl_intel_thread -lgfortran
#LIBS = -I${MKLROOT}/include -liomp5 -lpthread -lm -ldl -llapack -lblas
#LIBS =  -Wl,--start-group ${MKLROOT}/lib/ia32/libmkl_intel.a ${MKLROOT}/lib/ia32/libmkl_intel_thread.a ${MKLROOT}/lib/ia32/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
#LIBS =  -liomp5 -lpthread -lm -ldl
LIBS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
#LIBS =  -I${MKLROOT}/include   -Wl,--start-group  ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group  -lpthread -lm   -ldl
#LIBS = -I${MKLROOT}/include -lmkl_core  -liomp5 -lpthread -lm -ldl
endif
ifeq ($(UNAME), Darwin)
LIBS = -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
CC2 = icc -openmp
endif


#LIBS = -L${MKLROOT}/lib/intel64_lin/  -Wl,-rpath,${MKLROOT}/lib/intel64_lin/ -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -liomp5
#CC2 = icc -openmp

CFLAGS2 =   -std=c++11 -O3# -I${MKLROOT}/include
# CFLAGS2 = -g  -std=c++11

#LDFLAGS = -L${MKLROOT}/lib/intel64  -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm -fopenmp -llapack -lblas
#LDFLAGS =   -lmkl_intel_lp64 -lmkl_core  -ldl -lpthread -lm -fopenmp -llapack -lblas  -lmkl_intel_thread
#LDFLAG = -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl


 #icc  -O2 -DMKL_ILP64  -I$HOME/programs/papi/include ex_performance.c -o ex_performance -lm -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -liomp5 -ldl


.PHONY: default all clean

default: $(TARGET)
all: default

debug: CFLAGS2 = -g  -O0 -std=c++11 -Wall
debug: $(TARGET)

OBJECTS2 = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.h)


%.o: %.cpp $(HEADERS)
	$(CC2) $(CFLAGS2)  -c $< -o $@


.PRECIOUS: $(TARGET) $(OBJECTS2)


post:
	rm -r PostProc; mkdir PostProc

sync:
	rsync -r ./remote/*.cpp ./remote/*.h ./remote/*.in ./remote/Makefile  .

$(TARGET):    $(OBJECTS2)
	$(CC2)  $(OBJECTS2) -o $@   $(LIBS)

clean:
	-rm -f *.o
	-rm -f $(TARGET)
