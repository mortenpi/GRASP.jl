CC_FLAGS=-g
FC_FLAGS=-g -fPIC -std=f2008

TEST_BINDIR=../test/binaries

.PHONY: all
all: libgrasp.so $(TEST_BINDIR)/test-libgrasp

src/mixread.o: src/g2k_c_binding.o

libgrasp.so: src/mixread.o src/g2k_c_binding.o
	$(FC) -shared -o $@ $^

$(TEST_BINDIR)/test-libgrasp: test/libgrasp.o
	$(CXX) -o $@ $(CC_FLAGS) -Wl,-rpath,'$$ORIGIN/../../deps/' -L. $^ -lgrasp

%.o: %.f90
	$(FC) $(FC_FLAGS) -c -o $@ $^

%.o: %.cc
	$(CXX) $(CC_FLAGS) -c -o $@ $^

clean:
	-rm -vf libgrasp.so $(TEST_BINDIR)/test-libgrasp
	-rm -vf test/*.o src/*.o *.mod
