CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra
TARGET = ccalc.exe
SRCS = main.cpp bignum.cpp value.cpp parser.cpp evaluator.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp ccalc.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	del /Q $(OBJS) $(TARGET) 2>nul

.PHONY: all clean
