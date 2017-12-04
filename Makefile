FLAGS = -Wall
CFLAGS = -std=gnu11
CXXFLAGS = -std=gnu++11
CLIBS = -lm -lpng
BIN = ./bin
OBJ = ./obj
objects = pngUtil.o main.o
OBJS = $(objects:%.o=$(OBJ)/%.o)
binaries = main
BINS = $(binaries:%=$(BIN)/%)

ifeq (1, $(DEBUG))
FLAGS := $(FLAGS) -g
endif

.PHONY: all
all: $(BINS) 

$(BIN)/%: $(OBJ)/%.o $(OBJS)
	@echo Linking $@
	@mkdir -p $(BIN)
	g++ -o $@ $(OBJS) $(FLAGS) $(CXXFLAGS) $(CLIBS)

.PRECIOUS: $(OBJ)/%.o
$(OBJ)/%.o: ./src/%.c
	@echo Compiling $<
	@mkdir -p $(OBJ)
	gcc -MMD -c -o $@ $< $(FLAGS) $(CFLAGS)

$(OBJ)/%.o: ./src/%.cpp
	@echo Compiling $<
	@mkdir -p $(OBJ)
	g++ -MMD -c -o $@ $< $(FLAGS) $(CXXFLAGS)

.PHONY:clean
clean:
	@rm -f *.png
	@rm -f -r $(OBJ)
	@rm -f -r $(BIN)

-include $(OBJ)/*.d
