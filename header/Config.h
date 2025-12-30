#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

/**
 * Config class - Reads and stores configuration parameters from nmr.in file
 * 
 * File format:
 * Line 1: Input data filename
 * Line 2: Baseline adjustment
 * Line 3: Tolerance for numerical algorithms
 * Line 4: Filter type (0=none, 1=boxcar, 2=SG)
 * Line 5: Filter size (odd number)
 * Line 6: Number of filter passes
 * Line 7: Integration technique (0=Newton-Cotes, 1=Romberg, 2=Adaptive, 3=Quadrature)
 * Line 8: Output filename
 */
class Config {
public:
    string inputFilename;
    double baselineAdjustment;
    double tolerance;
    int filterType;  // 0=none, 1=boxcar, 2=SG
    int filterSize;
    int filterPasses;
    int integrationType;  // 0=Newton-Cotes, 1=Romberg, 2=Adaptive, 3=Quadrature
    string outputFilename;
    
    Config();
    bool readFromFile(const string& configFile);
    void print() const;
    
    string getFilterTypeName() const;
    string getIntegrationTypeName() const;
};

#endif // CONFIG_H
