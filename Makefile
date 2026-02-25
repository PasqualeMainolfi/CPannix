CC = clang
CFLAGS = -Wall -Wextra -g -O3

CPNX = cpannix/cpnx.c
HPNX = cpannix/cpnx.h
OPNX = target/temp/cpnx.o
OTEST = target/temp/pnx_test.o
CTEST = pnx_test.c

OBJ = $(OPNX) $(OTEST)
OUT = target/temp/cpnx_test
BUILD_FOLDER = target

ifeq ($(OS), Window_NT)
	MKDIR_TARGET = if not exist mkdir target mkdir target
	MKDIR_TEMP = if not exist mkdir target\\temp mkdir target\\temp
	RM_FILES = del /Q
	RM_DIR = del /S /Q
	EXE = .exe
	SHARED_EXT = dll
	SHARED_FLAG = -shared -BUILDING_CPNX
	SHARED_NAME = cpnx.$(SHARED_EXT)
else
	MKDIR_TARGET = mkdir -p target
	MKDIR_TEMP = mkdir -p target/temp
	RM_FILES = rm -f
	RM_DIR = rm -rf
	EXE =
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S), Darwin)
		SHARED_EXT = dylib
	else
		SHARED_EXT = so
	endif
	SHARED_FLAG = -fPIC -shared
	SHARED_NAME = libcpnx.$(SHARED_EXT)
endif

OUT := $(OUT)$(EXE)
SHARED_OUT = $(BUILD_FOLDER)/$(SHARED_NAME)

.PHONY: all clean clean_target clean_local test run_test shared
all: test

# --- BUILD SHARED LIB ---

shared: $(SHARED_OUT)

$(SHARED_OUT):
	$(MKDIR_TARGET)
	$(CC) $(CFLAGS) $(SHARED_FLAG) $(CPNX) -o $(SHARED_OUT)

# --- END BUILD SHARED LIB ---

# --- LOCAL TEST ---

test: $(OUT)

$(OUT): $(OBJ)
	$(CC) $(OBJ) -o $(OUT)

$(OPNX): $(CPNX) $(HPNX)
	$(MKDIR_TARGET)
	$(MKDIR_TEMP)
	$(CC) $(CFLAGS) -c $(CPNX) -o $(OPNX)

$(OTEST): $(CTEST) $(HPNX)
	$(MKDIR_TARGET)
	$(MKDIR_TEMP)
	$(CC) $(CFLAGS) -c $(CTEST) -o $(OTEST)

run_test: $(OUT)
	./$(OUT)

# --- END LOCAL TEST ---

clean_target:
	$(RM_DIR) target

clean_test:
	$(RM_FILES) target/temp/*.o target/temp/*.exe target/temp/cpnx_test target/*.dylib target/*.so target/*.dll

clean:
	$(RM_DIR) target
	$(RM_FILES) target/temp/*.o target/temp/*.exe target/temp/cpnx_test target/*.dylib target/*.so target/*.dll
