CXX = g++
CXXFLAGS = `root-config --cflags --libs`
TARGET = bin/main
SRC = main.cpp

all: $(TARGET)
	@./$(TARGET)

$(TARGET): $(SRC)
	@mkdir -p bin
	@$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

.PHONY: all
