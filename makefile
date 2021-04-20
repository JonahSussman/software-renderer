# 64 bit compilation target
CC = x86_64-w64-mingw32-g++
# 32 bit compilation target
# CC = i686-w64-mingw32

CFLAGS = -std=c++17 -O2
TFLAGS = -std=c++17 -O0 -ggdb -static-libstdc++ -static-libgcc
LFLAGS = -static-libstdc++ -static-libgcc -lmingw32 -lSDL2main -lSDL2

OBJS = src/main.cpp
TARGET = bin/software-renderer.exe

TOBJS = sw_math math_test

TTARGET = bin/software-renderer-tests.exe

ODIR = obj
SDIR = src

INCLUDE_PATHS = -Ilib/SDL2-2.0.14/x86_64-w64-mingw32/include/SDL2
LIBRARY_PATHS = -Llib/SDL2-2.0.14/x86_64-w64-mingw32/lib

all: $(OBJS)
	$(CC) $(OBJS) $(INCLUDE_PATHS) $(LIBRARY_PATHS) $(CFLAGS) $(LFLAGS) -o $(TARGET)

tests:
	$(CC) $(TFLAGS) src/math_test.cpp -o $(TTARGET)

# tests: $(TOBJS)
# 	$(CC) $(addsuffix .o, $(addprefix obj/, $(TOBJS))) $(INCLUDE_PATHS) $(LIBRARY_PATHS) $(TFLAGS) $(LFLAGS) -o $(TTARGET)

# $(TOBJS): 
# 	$(CC) -c $(addsuffix .cpp, $(addprefix src/, $@)) $(INCLUDE_PATHS) $(LIBRARY_PATHS) $(TFLAGS) $(LFLAGS) -o $(addsuffix .o, $(addprefix obj/, $@))