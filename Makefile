EXE_NAME := wave.exe

CXX := g++

SRC_EXTENSION := cpp

PLOG_DIR := ./external/plog
PLOG_INCLUDE := ${PLOG_DIR}/include/

PICOJSON_DIR := ./external/picojson
PICOJSON_INCLUDE := ${PICOJSON_DIR}

MUPARSER_VERSION := 2.2.5
MUPARSER_DIR := ./external/muparser-${MUPARSER_VERSION}
MUPARSER_INCLUDE := ${MUPARSER_DIR}/include/
MUPARSER_LIB := ${MUPARSER_DIR}/lib/

SUBDIRS := ${MUPARSER_DIR}
LIB_DIR := ${MUPARSER_LIB}

COMPILE_FLAGS = -Wall -Wextra -std=c++11 -O2 -g -Werror

INCLUDES := -I ${PLOG_INCLUDE} -I ${PICOJSON_INCLUDE} -I ${MUPARSER_INCLUDE}
DEFINES := -D RealT=double
LDFLAGS := -lrt -lmuparser

SHELL := /bin/bash

SRC_PATH = ./
BUILD_PATH := ./

# Will put the most recently modified sources first.
SOURCES := $(shell find ./	-maxdepth 1 -name '*.$(SRC_EXTENSION)' -printf '%T@\t%p\n' \
	 						| sort -k 1nr | cut -f2-)	

OBJECTS := $(SOURCES:%.$(SRC_EXTENSION)=%.$(SRC_EXTENSION).o)

.PHONY: subdirs ${SUBDIRS}

subdirs:
	$(MAKE) -C ${SUBDIRS}

$(EXE_NAME): $(OBJECTS) subdirs
	@echo "Linking: $@"
	$(CXX) $(OBJECTS) -o $@ -L ${LIB_DIR} $(LDFLAGS)

$(BUILD_PATH)/%.$(SRC_EXTENSION).o: %.$(SRC_EXTENSION)
	@echo $(SOURCES)
	@echo "Compiling $< -> $@"
	$(CXX) $(DEFINES) $(COMPILE_FLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -fv $(EXE_NAME) $(OBJECTS)
	rm -fv ${MUPARSER_DIR}/*.o ${MUPARSER_DIR}/lib/libmuparser*.so*
