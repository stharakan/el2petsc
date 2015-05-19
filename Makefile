include $(EL_DIR)/conf/ElVars
CPP_FLAGS = -openmp -O3 -std=c++11  

PSC_INC = -I$(PETSC_DIR)/include
PSC_LIB = -L$(PETSC_DIR)/lib -lpetsc


ALL_INCS = $(EL_LINK_FLAGS) $(PSC_INC) 
ALL_LIBS = -L./ $(EL_LIBS) $(PSC_LIB)

ELPSC_OBJ = el_petsc_utils.o
ELPSC_SRC = el_petsc_utils.cpp
ELPSC_DEPS = el_petsc_utils.hpp

TEST_BIN = test_funcs.exe
TEST_OBJ = test_funcs.o
TEST_SRC = test_funcs.cpp

ALL_OBJS = $(TEST_OBJ) $(ELPSC_OBJ)

all : el2petsc test

el2petsc : $(ELPSC_OBJ)

$(ELPSC_OBJ) : $(ELPSC_SRC) $(ELPSC_DEPS)
	$(CXX) $(EL_COMPILE_FLAGS) $(CPP_FLAGS) -c $(ELPSC_SRC) $(PSC_INC) $(EL_LINK_FLAGS) $(PSC_LIB) $(EL_LIBS) -o $@

test : $(TEST_BIN)

$(TEST_BIN) : $(ALL_OBJS)
	$(CXX) $(EL_COMPILE_FLAGS) $(CPP_FLAGS) $^ $(ALL_INCS) $(ALL_LIBS) -o $(TEST_BIN)

$(TEST_OBJ) : $(TEST_SRC) 
	$(CXX) $(EL_COMPILE_FLAGS) $(CPP_FLAGS) -c $(TEST_SRC) $(ALL_INCS) $(ALL_LIBS) -o $@

clean:
	rm *.o
	rm *.exe

