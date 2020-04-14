TARGET = main
UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
CC2 = icc 
#LIBS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
LIBS = -lm -ldl
endif
ifeq ($(UNAME), Darwin)
#LIBS = -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
LIBS = -lm -ldl
CC2 = icc 
endif

CFLAGS2 = -std=c++11 -O3

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
