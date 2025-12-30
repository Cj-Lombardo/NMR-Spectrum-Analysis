# Makefile for NMR Spectrum Analysis
# CSC/MAT 335 - Project 1
# Place this file in the build/ directory

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2 -I../header
LDFLAGS = -larmadillo

# Directories
SRC_DIR = ../src
HEADER_DIR = ../header
DATA_DIR = ../data

# Target executable
TARGET = nmr_analysis

# Object files (in build directory)
OBJECTS = main.o Config.o DataReader.o Filter.o CubicSpline.o Integration.o PeakDetector.o DataWriter.o

# Header files (for dependencies)
HEADERS = $(HEADER_DIR)/Config.h \
          $(HEADER_DIR)/DataReader.h \
          $(HEADER_DIR)/Filter.h \
          $(HEADER_DIR)/CubicSpline.h \
          $(HEADER_DIR)/Integration.h \
          $(HEADER_DIR)/PeakDetector.h \
          $(HEADER_DIR)/DataWriter.h

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJECTS)
	@echo "Linking $(TARGET)..."
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) $(LDFLAGS)
	@echo "Build complete!"

# Compile main.cpp
main.o: $(SRC_DIR)/main.cpp $(HEADERS)
	@echo "Compiling main.cpp..."
	$(CXX) $(CXXFLAGS) -c $(SRC_DIR)/main.cpp -o main.o

# Compile Config.cpp
Config.o: $(SRC_DIR)/Config.cpp $(HEADER_DIR)/Config.h
	@echo "Compiling Config.cpp..."
	$(CXX) $(CXXFLAGS) -c $(SRC_DIR)/Config.cpp -o Config.o

# Compile DataReader.cpp
DataReader.o: $(SRC_DIR)/DataReader.cpp $(HEADER_DIR)/DataReader.h
	@echo "Compiling DataReader.cpp..."
	$(CXX) $(CXXFLAGS) -c $(SRC_DIR)/DataReader.cpp -o DataReader.o

# Compile Filter.cpp
Filter.o: $(SRC_DIR)/Filter.cpp $(HEADER_DIR)/Filter.h
	@echo "Compiling Filter.cpp..."
	$(CXX) $(CXXFLAGS) -c $(SRC_DIR)/Filter.cpp -o Filter.o

# Compile CubicSpline.cpp
CubicSpline.o: $(SRC_DIR)/CubicSpline.cpp $(HEADER_DIR)/CubicSpline.h
	@echo "Compiling CubicSpline.cpp..."
	$(CXX) $(CXXFLAGS) -c $(SRC_DIR)/CubicSpline.cpp -o CubicSpline.o

# Compile Integration.cpp
Integration.o: $(SRC_DIR)/Integration.cpp $(HEADER_DIR)/Integration.h $(HEADER_DIR)/CubicSpline.h
	@echo "Compiling Integration.cpp..."
	$(CXX) $(CXXFLAGS) -c $(SRC_DIR)/Integration.cpp -o Integration.o

# Compile PeakDetector.cpp
PeakDetector.o: $(SRC_DIR)/PeakDetector.cpp $(HEADER_DIR)/PeakDetector.h $(HEADER_DIR)/CubicSpline.h $(HEADER_DIR)/Integration.h
	@echo "Compiling PeakDetector.cpp..."
	$(CXX) $(CXXFLAGS) -c $(SRC_DIR)/PeakDetector.cpp -o PeakDetector.o

# Compile DataWriter.cpp
DataWriter.o: $(SRC_DIR)/DataWriter.cpp $(HEADER_DIR)/DataWriter.h $(HEADER_DIR)/CubicSpline.h $(HEADER_DIR)/PeakDetector.h
	@echo "Compiling DataWriter.cpp..."
	$(CXX) $(CXXFLAGS) -c $(SRC_DIR)/DataWriter.cpp -o DataWriter.o

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(OBJECTS) $(TARGET) *.txt
	@echo "Clean complete!"

# Run the program with default config from data directory
run: $(TARGET)
	./$(TARGET) $(DATA_DIR)/nmr.in

# Help target
help:
	@echo "NMR Spectrum Analysis - Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all       - Build the program (default)"
	@echo "  clean     - Remove build artifacts"
	@echo "  run       - Build and run with config from data/"
	@echo "  help      - Show this help message"

.PHONY: all clean run help
