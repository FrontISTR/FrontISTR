#> ggtools Makefile

##> compiler setting
FC     = mpif90
FFLAGS = -fPIC -O2 -mtune=native -march=native -std=legacy -Wno-missing-include-dirs
CC     = mpicc
CFLAGS = -fPIC -O2

##> directory setting
MOD_DIR = -J ./include
INCLUDE = -I /usr/include -I ./include -I ./submodule/monolis_utils/include
USE_LIB = -L./lib -lggtools -L./submodule/monolis_utils/lib -lmonolis_utils
BIN_DIR = ./bin
SRC_DIR = ./src
TST_DIR = ./src_test
OBJ_DIR = ./obj
LIB_DIR = ./lib
WRAP_DIR= ./wrapper
TST_WRAP_DIR = ./wrapper_test
DRV_DIR = ./driver
LIBRARY = libggtools.a
CPP     = -cpp $(FLAG_DEBUG)

##> option setting
ifdef FLAGS
	comma:= ,
	empty:=
	space:= $(empty) $(empty)
	DFLAGS = $(subst $(comma), $(space), $(FLAGS))

	ifeq ($(findstring DEBUG, $(DFLAGS)), DEBUG)
		FFLAGS  = -fPIC -O2 -std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow -Wno-missing-include-dirs -Wall
		CFLAGS  = -fPIC -O2 -g -ggdb
	endif

	ifeq ($(findstring INTEL, $(DFLAGS)), INTEL)
		FC      = mpiifort
		FFLAGS  = -fPIC -O2 -align array64byte  -nofor-main
		CC      = mpiicc
		CFLAGS  = -fPIC -O2 -no-multibyte-chars
		MOD_DIR = -module ./include
	endif

	ifeq ($(findstring SUBMODULE, $(DFLAGS)), SUBMODULE)
		INCLUDE = -I /usr/include -I ./include -I ../monolis_utils/include -I ../../include
		USE_LIB = -L./lib -lggtools -L../monolis_utils/lib -lmonolis_utils
	endif
endif

##> other commands
MAKE = make
CD   = cd
CP   = cp
RM   = rm -rf
AR   = - ar ruv

##> **********
##> target (1)
LIB_TARGET = $(LIB_DIR)/$(LIBRARY)

##> source file define
SRC_DEF = \
  def_element.f90 \
  def_global_bucket.f90 \
  def_local_bucket.f90 \
  distance_determination.f90

SRC_BUC = \
  element.f90 \
  bucket_util.f90 \
  global_bucket.f90 \
  local_bucket.f90

##> C wrapper section
#SRC_GRAPH_C =

SRC_ALL_C = \
$(addprefix graph/, $(SRC_GRAPH_C))

##> all targes
SRC_ALL = \
$(addprefix define/, $(SRC_DEF)) \
$(addprefix bucket/, $(SRC_BUC))

##> lib objs
LIB_SOURCES = \
$(addprefix $(SRC_DIR)/,  $(SRC_ALL)) \
$(addprefix $(WRAP_DIR)/, $(SRC_ALL_C)) \
./src/ggtools.f90
LIB_OBJSt   = $(subst $(SRC_DIR), $(OBJ_DIR), $(LIB_SOURCES:.f90=.o))
LIB_OBJS    = $(subst $(WRAP_DIR), $(OBJ_DIR), $(LIB_OBJSt:.c=.o))

##> **********
##> target (2) test for fotran
TEST_TARGET = $(TST_DIR)/ggtools_test

##> lib objs
TST_SRC_ALL = $(SRC_ALL) ggtools.f90
TST_SOURCES = $(addprefix $(TST_DIR)/, $(TST_SRC_ALL))
TST_OBJSt   = $(subst $(TST_DIR), $(OBJ_DIR), $(TST_SOURCES:.f90=_test.o))
TST_OBJS    = $(TST_OBJSt:.c=_test.o)

##> **********
##> target (3) test for fotran
TEST_C_TARGET = $(TST_WRAP_DIR)/ggtools_c_test

##> lib objs
#SRC_GRAPH_C_TEST =

SRC_ALL_C_TEST = \
$(addprefix graph/, $(SRC_GRAPH_C_TEST))

TST_SRC_C_ALL = $(SRC_ALL_C_TEST) ggtools_c_test.c
TST_C_SOURCES = $(addprefix $(TST_WRAP_DIR)/, $(TST_SRC_C_ALL))
TST_C_OBJS    = $(subst $(TST_WRAP_DIR), $(OBJ_DIR), $(TST_C_SOURCES:.c=.o))

##> target
all: \
	cp_header \
	$(LIB_TARGET) \
	$(TEST_TARGET) \
	$(TEST_C_TARGET)

lib: \
	cp_header \
	$(LIB_TARGET)

$(LIB_TARGET): $(LIB_OBJS)
	$(AR) $@ $(LIB_OBJS)

$(TEST_TARGET): $(TST_OBJS)
	$(FC) $(FFLAGS) -o $@ $(TST_OBJS) $(USE_LIB)

$(TEST_C_TARGET): $(TST_C_OBJS)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(TST_C_OBJS) $(USE_LIB)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

$(OBJ_DIR)/%.o: $(TST_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

$(OBJ_DIR)/%.o: $(WRAP_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

$(OBJ_DIR)/%.o: $(WRAP_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

$(OBJ_DIR)/%.o: $(TST_WRAP_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

cp_header:
#	$(CP) ./wrapper/ggtools.h ./include/

clean:
	$(RM) \
	$(LIB_OBJS) \
	$(TST_OBJS) \
	$(TST_C_OBJS) \
	$(LIB_TARGET) \
	$(TEST_TARGET) \
	$(TEST_C_TARGET) \
	./include/*.h \
	./include/*.mod \
	./bin/*

.PHONY: clean
