HOSTNAME := $(shell hostname)

# For reproducability compile the git key (hash) into the code
GITKEY = $(shell git log | head -1 | sed -e 's/commit //g')

EXECUTABLE = impera
all: ${EXECUTABLE}

OBJ = \
	control.o \
	warnings.o \
	impera.o \
	window.o \
	control.o \
	main.o

ifdef ENABLEGPU
	OBJ += GPUdrivers.o
endif


.PHONY: clean allclean
clean:
	-rm -f *.o *~ *.d *.bmp

allclean: clean
	-rm -f ${EXECUTABLE}

# CXX    = icpc
CXX    = g++

SRC ?= ./src
INC ?=
# INC += -I/usr/local/cuda/include
INC += -I../include

LDFLAGS ?=
ERRFLAGS ?=
FEAFLAGS ?=

FEAFLAGS += -D_GIT_KEY=$(GITKEY) 
# ERRFLAGS = -Werror -Wall -Wfatal-errors
# FEAFLAGS += -fopenmp -pthread
# FEAFLAGS += -fopenmp
# FEAFLAGS += -fno-omit-frame-pointer
# FEAFLAGS += -D NO_UNIT_TESTS ## for production

## include also code which is not meant for release,
## use github.com/BR903/cppp -UDEVEL to do partial preprocessing
FEAFLAGS += -D DEVEL 
FEAFLAGS += -D HAS_WINDOW -D GL_SILENCE_DEPRECATION
LDFLAGS += -framework GLUT -framework OpenGL

# normal mode
OPTFLAGS = -Ofast
# OPTFLAGS = -O0
### -mtune=power8 -mcpu=power8 -mpower8-fusion -mpower8-vector -mvsx -maltivec -mdirect-move -fstrict-aliasing -fomit-frame-pointer -fno-schedule-insns
# OPTFLAGS = -O2
OPTFLAGS += -g -pedantic -Wall -Warray-bounds -Wno-sign-compare -Wno-char-subscripts -Wno-format
# OPTFLAGS += -Wno-unused-variable
# OPTFLAGS += -Wno-nested-anon-types


CUOPTFLAGS =
#### production mode
CUOPTFLAGS += -O3 --use_fast_math
#### debug mode
# OPTFLAGS  = -O0 -g
# CUOPTFLAGS += -g -O0 -G
# CUOPTFLAGS += -DDEBUGGPU

EXTFLAGS ?=

ifeq ($(HOSTNAME),MacBook-Pro.local)
# ifeq ($(HOSTNAME),zam044)
# 	EXTFLAGS += -D HAS_no_MKL
	LDFLAGS += -lm -lblas -llapack
else
	## eq ($(HOSTNAME),zam054)
	### works with g++ on zam054 (JSC workstation)
	LDFLAGS += -lm -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread 
endif

FLAGS = ${ERRFLAGS} ${FEAFLAGS} ${OPTFLAGS} ${EXTFLAGS}

CCFLAGS  = -std=gnu99 ${FLAGS}
CXXFLAGS = -std=c++11 ${FLAGS}

DEP=$(OBJ:%.o=%.d)

${EXECUTABLE}: ${OBJ}
	${CXX} ${CXXFLAGS} $^ ${LDFLAGS}  -o $@

%.o:${SRC}/%.c %.d
	${CC} ${CCFLAGS} ${INC} -c $<
%.o:${SRC}/%.cxx %.d
	${CXX} ${CXXFLAGS} ${INC} -c $<

%.o:%.c %.d
	${CC} ${CCFLAGS} ${INC} -c $<
%.o:%.cxx %.d
	${CXX} ${CXXFLAGS} ${INC} -c $<

%.d:${SRC}/%.c
	${CC} ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:${SRC}/%.cxx
	${CXX} ${CXXFLAGS} ${INC} -MM $^ -MF $@

%.d:%.c
	${CC} ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:%.cxx
	${CXX} ${CXXFLAGS} ${INC} -MM $^ -MF $@

-include ${DEP}
