all: test-out

test: test.cpp matrix.h util.h matrix_traits.h
	reset && clang++-3.9 -std=c++1z -o test test.cpp -O3 
test-out: test FORCE
	./test | tee test-out

clean:
	rm -rf ./test ./test-out

FORCE: ;
