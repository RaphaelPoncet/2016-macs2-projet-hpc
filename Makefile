EXE_NAME := wave.exe

CXX := g++

SRC_EXTENSION := cpp

COMPILE_FLAGS = -Wall -Wextra -Werror -std=c++11 -O2 -g

INCLUDES := -I ./external/plog/include -I ./external/picojson/
DEFINES := -D RealT=double
LDFLAGS := -lrt


SHELL := /bin/bash

SRC_PATH = ./
BUILD_PATH := ./

# Will put the most recently modified sources first.
SOURCES := $(shell find ./	-maxdepth 1 -name '*.$(SRC_EXTENSION)' -printf '%T@\t%p\n' \
	 						| sort -k 1nr | cut -f2-)	

OBJECTS := $(SOURCES:%.$(SRC_EXTENSION)=%.$(SRC_EXTENSION).o)

$(EXE_NAME): $(OBJECTS) 
	@echo "Linking: $@"
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

$(BUILD_PATH)/%.$(SRC_EXTENSION).o: %.$(SRC_EXTENSION)
	@echo $(SOURCES)
	@echo "Compiling $< -> $@"
	$(CXX) $(DEFINES) $(COMPILE_FLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(EXE_NAME) $(OBJECTS)
