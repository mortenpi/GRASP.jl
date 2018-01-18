CC_FLAGS=-g
FC_FLAGS=-g -fPIC -std=f2008

.PHONY: all
all: libgrasp.so ../test/test-libgrasp

libgrasp.so: mixread.o
	$(FC) -shared -o $@ $^

../test/test-libgrasp: test-libgrasp.o
	$(CXX) -o $@ $(CC_FLAGS) -Wl,-rpath,'$$ORIGIN/../deps/' -L. $^ -lgrasp

%.o: %.f90
	$(FC) $(FC_FLAGS) -c -o $@ $^

%.o: %.cc
	$(CXX) $(CC_FLAGS) -c -o $@ $^

clean:
	-rm -vf libgrasp.so ../test/test-libgrasp
	-rm -vf *.o *.mod
