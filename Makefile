
OBJECT_FILES = gen.o coo.o grid.o physics.o openblaswrapper.o lapackwrapper.o rsbwrapper.o precond.o splrep.o splopb4dcg.o
HEADER_FILES = ./src/gen.h ./src/coo.h ./src/grid.h ./src/physics.h ./src/openblaswrapper.h ./src/lapackwrapper.h ./src/rsbwrapper.h ./src/precond.h ./src/splrep.h ./src/splopb4dcg.h
GEN_PATH = include/gen.h
OPENBLAS_PATH1 = $HOME/OpenBlas/lib
OPENBLAS_PATH2 = $HOME/OpenBlas/include
LAPACK_PATH = $HOME/ARPACK/LAPACK
RSB_PATH = $HOME/librsb/lib

I_OPTS = "{includedir}"
includedir = "${prefix}/include"
prefix = "/nfs/users3/philst/librsb"

CFLAGS = -std=c11 `librsb-config --I_opts --cflags` -O2 -Wall -Wextra -Isrc -rdynamic -pedantic $(OPTFLAGS) -march=native -mfpmath=sse -msse2 -fopenmp
MATHFLAGS = -lm
#LIBS = -ldl $(OPTLIBS)
#PREFIX ?= /usr/local

#SOURCES = $(wildcard src/**/*.c src/*.c)
#OBJECTS = $(patsubst %.c,%.o,%(SOURCES))

#TARGET = build/libLREP.a
#SO_TARGET = $(patsubst %.a,%.so,$(TARGET))

# The Target Build
#all: $(TARGET) $(SO_TARGET)

#dev: CFLAGS=-g -WALL -Isrc -Wall -Wextra $(OPTFLAGS)
#dev: all

#$(TARGET): CFLAGS += -fPIC
#$(TARGET): build $(OBJECTS)
#	ar rcs $@ $(OBJECTS)
#	ranlib $@

#$(SO_TARGET): $(TARGET) $(OBJECTS)
#	$(CC) -shared -o $@ $(OBJECTS)

#build:
#	@mkdir -p build
#	@mkdir -p bin

# The Cleaner
#clean:
#	rm -rf build $(OBJECTS) $(TESTS)
#	rm -f test/tests.log
#	find . -name "*.gc*" -exec rm {} \;
#	rm -rf 'find . -name "*.dSYM" -print'

all: liblrep.a

liblrep.a: $(OBJECT_FILES) $(HEADER_FILES)
	ar r liblrep.a $(OBJECT_FILES)
	ranlib liblrep.a

gen.o: ./src/gen.c ./src/gen.h
	gcc $(CFLAGS) -c ./src/gen.c

coo.o: ./src/coo.c ./src/coo.h
	gcc $(CFLAGS) -c ./src/coo.c

grid.o: ./src/grid.c ./src/grid.h
	gcc $(CFLAGS) -c ./src/grid.c

physics.o: ./src/physics.c ./src/physics.h
	gcc $(CFLAGS) -c ./src/physics.c

openblaswrapper.o: ./src/openblaswrapper.c ./src/openblaswrapper.h
	gcc $(CFLAGS) -c ./src/openblaswrapper.c

lapackwrapper.o: ./src/lapackwrapper.c ./src/lapackwrapper.h
	gcc $(CFLAGS) -c ./src/lapackwrapper.c

rsbwrapper.o: ./src/rsbwrapper.c ./src/rsbwrapper.h
	gcc $(CFLAGS) -c ./src/rsbwrapper.c

precond.o: ./src/precond.c ./src/precond.h
	gcc $(CFLAGS) -c ./src/precond.c

splrep.o: ./src/splrep.c ./src/splrep.h
	gcc $(CFLAGS) -c ./src/splrep.c

splopb4dcg.o: ./src/splopb4dcg.c ./src/splopb4dcg.h
	gcc $(CFLAGS) -c ./src/splopb4dcg.c

clean:
	rm *.a *.o *.ex

test_diff2.ex: ./test_rsb_setup_diff2.c
	gcc $(CFLAGS) -march=native -mfpmath=sse -msse2 -c test_rsb_setup_diff2.c
	gcc -std=c11 -g test_rsb_setup_diff2.o -o test.ex -L. -llrep `librsb-config --static --ldflags --extra_libs` -L$(LAPACK_PATH) -llapack -I$(OPENBLAS_PATH2) -L$(OPENBLAS_PATH1) -lopenblas -pg -lm

test_fourier2.ex: ./test_fourier2.c
	gcc $(CFLAGS) -march=native -mfpmath=sse -msse2 -c test_fourier2.c
	gcc -std=c11 -g test_fourier2.o -o test.ex -L. -llrep `librsb-config --static --ldflags --extra_libs` -L$(LAPACK_PATH) -llapack -I$(OPENBLAS_PATH2) -L$(OPENBLAS_PATH1) -lopenblas -pg -lm

test_kron.ex: ./test_kron.c
	gcc $(CFLAGS) -march=native -mfpmath=sse -msse2 -c test_kron.c
	gcc -std=c11 -g test_kron.o -o test.ex -L. -llrep `librsb-config --static --ldflags --extra_libs` -L$(LAPACK_PATH) -llapack -I$(OPENBLAS_PATH2) -L$(OPENBLAS_PATH1) -lopenblas -pg -lm

test_grid.ex: ./test_grid.c
	gcc $(CFLAGS) -march=native -mfpmath=sse -msse2 -c ./test_grid.c
	gcc -std=c11 -g test_grid.o -o test.ex -L. -llrep `librsb-config --static --ldflags --extra_libs` -L$(LAPACK_PATH) -llapack -I$(OPENBLAS_PATH2) -L$(OPENBLAS_PATH1) -lopenblas -pg -lm

test_trap.ex: ./test_trap.c
	gcc $(CFLAGS) -march=native -mfpmath=sse -msse2 -c ./test_trap.c
	gcc -std=c11 -g test_trap.o -o test.ex -L. -llrep `librsb-config --static --ldflags --extra_libs` -L$(LAPACK_PATH) -llapack -I$(OPENBLAS_PATH2) -L$(OPENBLAS_PATH1) -lopenblas -pg -lm

test_uint.ex: ./test_uint.c
	gcc $(CFLAGS) -march=native -mfpmath=sse -msse2 -c ./test_uint.c
	gcc -std=c11 -g test_uint.o -o test.ex -L. -llrep `librsb-config --static --ldflags --extra_libs` -L$(LAPACK_PATH) -llapack -I$(OPENBLAS_PATH2) -L$(OPENBLAS_PATH1) -lopenblas -pg -lm


valgrind: 
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind.out.txt ./main.ex

