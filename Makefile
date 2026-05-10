CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

TARGET = ccalc.exe
GRAPH_TARGET = ccalc_graph.exe

SRCS = main.cpp bignum.cpp value.cpp parser.cpp evaluator.cpp
OBJS = $(SRCS:.cpp=.o)

IMGUI_DIR = imgui-1.92.7
IMGUI_SRCS = $(IMGUI_DIR)/imgui.cpp \
             $(IMGUI_DIR)/imgui_draw.cpp \
             $(IMGUI_DIR)/imgui_tables.cpp \
             $(IMGUI_DIR)/imgui_widgets.cpp \
             $(IMGUI_DIR)/imgui_demo.cpp \
             $(IMGUI_DIR)/backends/imgui_impl_win32.cpp \
             $(IMGUI_DIR)/backends/imgui_impl_dx11.cpp
IMGUI_OBJS = $(IMGUI_SRCS:.cpp=.o)

CALC_SRCS = bignum.cpp value.cpp parser.cpp evaluator.cpp
CALC_OBJS = $(CALC_SRCS:.cpp=.o)

GRAPH_SRCS = graph.cpp
GRAPH_OBJS = $(GRAPH_SRCS:.cpp=.o)

IMGUI_CXXFLAGS = -I$(IMGUI_DIR) -I$(IMGUI_DIR)/backends

LIBS_GRAPH = -ld3d11 -ld3dcompiler -ldxgi -ldwmapi -lgdi32 -luser32 -lshell32

all: $(TARGET) $(GRAPH_TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp ccalc.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(IMGUI_DIR)/%.o: $(IMGUI_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(IMGUI_CXXFLAGS) -c -o $@ $<

$(GRAPH_TARGET): $(GRAPH_OBJS) $(CALC_OBJS) $(IMGUI_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS_GRAPH)

graph.o: graph.cpp ccalc.h
	$(CXX) $(CXXFLAGS) $(IMGUI_CXXFLAGS) -c -o $@ $<

clean:
	del /Q $(OBJS) $(GRAPH_OBJS) $(CALC_OBJS) $(IMGUI_OBJS) $(TARGET) $(GRAPH_TARGET) 2>nul

.PHONY: all clean
