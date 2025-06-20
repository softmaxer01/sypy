CXX = g++
CXXFLAGS = -std=c++11 -Iinclude
SRC = src/main.cpp src/sypy.cpp
TARGET = main

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET) 