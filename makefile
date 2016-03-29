CC=/usr/bin/g++
# CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11 -march=native -pipe
LDFLAGS= -pthread -pipe -flto


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O3  -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O3 -g
LDFLAGS=-g
endif

ifeq ($(polly),1)
CFLAGS+= -Xclang -load -Xclang /local/polly/llvm_build/lib/LLVMPolly.so -mllvm -polly -mllvm -polly-vectorizer=stripmine
CC=/local/polly/llvm_build/bin/clang
endif


# Directory organisation
BASEDIR = .
SRC = $(BASEDIR)/src
BUILD = $(BASEDIR)/build
BIN = $(BASEDIR)/bin


# Program name
TARGET = ftl


# Objects names
OBJECTS = $(BUILD)/main.o $(BUILD)/distances.o $(BUILD)/mapping.o $(BUILD)/ssw_cpp.o $(BUILD)/ssw.o $(BUILD)/utils.o


all: init $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(BIN)/$(TARGET) $^

$(BUILD)/%.o: $(SRC)/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

$(BUILD)/%.o: $(SRC)/%.c
	$(CC) $(CFLAGS) -c -o $@ $^

clean:
	rm -rf $(BUILD)/*.o
	rm -rf $(BIN)/$(TARGET)

init: 
	mkdir -p $(BUILD) $(BIN)

rebuild: clean $(TARGET)
