CC_FLAGS=-g
FC_FLAGS=-g -fPIC -std=f2008

.PHONY: all
all: libgrasp.so ../test/test-libgrasp

libgrasp.so: mixread.o
	gfortran -shared -o $@ $^

../test/test-libgrasp: test-libgrasp.o
	g++ -o $@ $(CC_FLAGS) -Wl,-rpath,'$$ORIGIN/../deps/' -L. $^ -lgrasp

%.o: %.f90
	gfortran $(FC_FLAGS) -c -o $@ $^

%.o: %.cc
	g++ $(CC_FLAGS) -c -o $@ $^

clean:
	-rm -vf libgrasp.so ../test/test-libgrasp
	-rm -vf *.o *.mod
