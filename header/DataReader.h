#ifndef DATAREADER_H
#define DATAREADER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

/**
 * DataReader class - Reads NMR spectrum data from file
 * 
 * Data format: Two columns (x-value, y-value) separated by whitespace
 */
class DataReader {
public:
    vector<double> xData;
    vector<double> yData;
    
    DataReader();
    bool readFromFile(const string& filename);
    size_t size() const { return xData.size(); }
    void print(size_t numPoints = 10) const;
    
    // Sort data by x-values in ascending order (keeping x,y pairs together)
    void sortData();
    
    // Check if data is sorted in ascending order
    bool isSorted() const;
    
    // Find and shift data so TMS peak is at x=0.0
    double findAndShiftTMS(double baselineAdjustment);
    
    // Correct baseline by subtracting estimated baseline value
    // Returns the baseline value that was subtracted
    double correctBaseline();
};

#endif // DATAREADER_H
