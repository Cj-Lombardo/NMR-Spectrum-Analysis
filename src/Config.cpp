#include "Config.h"

using namespace std;

/**
 * Constructor - Initialize with default values
 */
Config::Config() 
    : inputFilename(""), baselineAdjustment(0.0), tolerance(1e-8),
      filterType(0), filterSize(0), filterPasses(0), integrationType(0),
      outputFilename("analysis.txt") {
}

/**
 * Read configuration from file
 * @param configFile - path to configuration file (typically nmr.in)
 * @return true if successful, false otherwise
 */
bool Config::readFromFile(const string& configFile) {
    ifstream inFile(configFile);
    
    if (!inFile.is_open()) {
        cerr << "Error: Cannot open configuration file: " << configFile << endl;
        return false;
    }
    
    string line;
    int lineNum = 0;
    
    // Read input filename (line 1)
    if (getline(inFile, line)) {
        istringstream iss(line);
        iss >> inputFilename;
        lineNum++;
    }
    
    // Read baseline adjustment (line 2)
    if (getline(inFile, line)) {
        istringstream iss(line);
        iss >> baselineAdjustment;
        lineNum++;
    }
    
    // Read tolerance (line 3)
    if (getline(inFile, line)) {
        istringstream iss(line);
        iss >> tolerance;
        lineNum++;
    }
    
    // Read filter type (line 4)
    if (getline(inFile, line)) {
        istringstream iss(line);
        iss >> filterType;
        lineNum++;
    }
    
    // Read filter size (line 5)
    if (getline(inFile, line)) {
        istringstream iss(line);
        iss >> filterSize;
        lineNum++;
    }
    
    // Read filter passes (line 6)
    if (getline(inFile, line)) {
        istringstream iss(line);
        iss >> filterPasses;
        lineNum++;
    }
    
    // Read integration type (line 7)
    if (getline(inFile, line)) {
        istringstream iss(line);
        iss >> integrationType;
        lineNum++;
    }
    
    // Read output filename (line 8)
    if (getline(inFile, line)) {
        istringstream iss(line);
        iss >> outputFilename;
        lineNum++;
    }
    
    inFile.close();
    
    if (lineNum < 8) {
        cerr << "Error: Configuration file incomplete. Expected 8 lines, found " 
                  << lineNum << endl;
        return false;
    }
    
    // Validate filter size is odd (if filtering is enabled)
    if (filterType != 0 && filterSize % 2 == 0) {
        cerr << "Warning: Filter size should be odd. Adjusting from " 
                  << filterSize << " to " << (filterSize + 1) << endl;
        filterSize++;
    }
    
    return true;
}

/**
 * Print configuration to console
 */
void Config::print() const {
    cout << "\n-=> NMR ANALYSIS <=-\n" << endl;
    cout << "Program Options" << endl;
    cout << "===============================" << endl;
    cout << "Input File          : " << inputFilename << endl;
    cout << "Baseline Adjustment : " << baselineAdjustment << endl;
    cout << "Tolerance           : " << tolerance << endl;
    cout << "Filter Type         : " << getFilterTypeName() << endl;
    if (filterType != 0) {
        cout << "Filter Size         : " << filterSize << endl;
        cout << "Filter Passes       : " << filterPasses << endl;
    }
    cout << "Integration Method  : " << getIntegrationTypeName() << endl;
    cout << "Output File         : " << outputFilename << endl;
    cout << endl;
}

/**
 * Get filter type name
 */
string Config::getFilterTypeName() const {
    switch (filterType) {
        case 0: return "None, Filtering is Off";
        case 1: return "Boxcar (Cyclic)";
        case 2: return "Savitzky-Golay";
        default: return "Unknown";
    }
}

/**
 * Get integration type name
 */
string Config::getIntegrationTypeName() const {
    switch (integrationType) {
        case 0: return "Newton-Cotes";
        case 1: return "Romberg";
        case 2: return "Adaptive Quadrature";
        case 3: return "Gauss-Legendre Quadrature";
        default: return "Unknown";
    }
}
