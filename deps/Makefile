CC_FLAGS=-g
FC_FLAGS=-g -fPIC -std=f2008 -cpp $(FC_FLAGS_USER)

TEST_BINDIR=../test/binaries

.PHONY: all
all: libgrasp.so \
	$(TEST_BINDIR)/test-libgrasp \
	$(TEST_BINDIR)/test-readrwfn \
	build_options.log

.PHONY: build_options.log
build_options.log:
	-echo "Build options at `date`:" >> $@
	-echo "  FC:       ${FC}" >> $@
	-echo "  CC_FLAGS: ${CC_FLAGS}" >> $@
	-echo "  FC_FLAGS: ${FC_FLAGS}" >> $@
	-echo >> $@

src/mixread.o: src/g2k_c_binding.o
src/rwfnread.o: src/g2k_c_binding.o

libgrasp.so: src/mixread.o src/rwfnread.o src/g2k_c_binding.o
	$(FC) -shared -o $@ $^

# Targets for test binaries.
# Note: they have an implicit dependence on libgrasp.so
$(TEST_BINDIR)/test-%: test/%.o
	$(CXX) -o $@ $(CC_FLAGS) -Wl,-rpath,'$$ORIGIN/../../deps/' -L. $^ -lgrasp

%.o: %.f90
	$(FC) $(FC_FLAGS) -c -o $@ $^

%.o: %.cc
	$(CXX) $(CC_FLAGS) -c -o $@ $^

clean:
	-rm -vf libgrasp.so $(TEST_BINDIR)/test-libgrasp
	-rm -vf test/*.o src/*.o *.mod
	-rm -vf build_options.log
