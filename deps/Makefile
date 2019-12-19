CC_FLAGS=-g
FC_FLAGS=-g -fPIC -std=f2008 -cpp $(FC_FLAGS_USER)

TEST_BINDIR=../test/binaries

.PHONY: all
all: libgrasp.so \
	libgraspjl.so \
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

libgraspjl.so: src/mixread.o src/rwfnread.o src/g2k_c_binding.o
	$(FC) -shared -o $@ $^

grasp:
	./getgrasp.sh

libgrasp.so: grasp
	cp grasp/lib/libgrasp.so libgrasp.so

# Targets for test binaries.
# Note: they have an implicit dependence on libgrasp.so
$(TEST_BINDIR)/test-%: test/%.o
	$(CXX) -o $@ $(CC_FLAGS) -Wl,-rpath,'$$ORIGIN/../../deps/' -L. $^ -lgraspjl

%.o: %.f90
	$(FC) $(FC_FLAGS) -c -o $@ $^

%.o: %.cc
	$(CXX) $(CC_FLAGS) -c -o $@ $^

clean:
	-rm -vf libgrasp.so
	-rm -vf libgraspjl.so $(TEST_BINDIR)/test-libgrasp
	-rm -vf grasp-path.jl
	-rm -vf test/*.o src/*.o *.mod
	-rm -vf *.log

clean-grasp:
	-rm -vf grasp/
