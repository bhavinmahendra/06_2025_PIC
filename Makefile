CXX = g++
CXXFLAGS = `root-config --cflags --libs`
LIB = /opt/homebrew/Cellar/boost/1.88.0/include/
TARGET = bin/main
SRC = main.cpp

all: $(TARGET)
	@./$(TARGET)

$(TARGET): $(SRC)
	@mkdir -p bin
	@$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) -I $(LIB)

.PHONY: all
