#include "DataReader.h"
#include <algorithm>
#include <cmath>

using namespace std;

/**
 * Constructor
 */
DataReader::DataReader() {
}

/**
 * Read data from file
 * @param filename - path to data file
 * @return true if successful, false otherwise
 */
bool DataReader::readFromFile(const string& filename) {
    ifstream inFile(filename);
    
    if (!inFile.is_open()) {
        cerr << "Error: Cannot open data file: " << filename << endl;
        return false;
    }
    
    xData.clear();
    yData.clear();
    
    string line;
    double x, y;
    
    while (getline(inFile, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        istringstream iss(line);
        if (iss >> x >> y) {
            xData.push_back(x);
            yData.push_back(y);
        }
    }
    
    inFile.close();
    
    cout << "Read " << xData.size() << " data points from " << filename << endl;
    
    if (xData.empty()) {
        cerr << "Error: No data points read from file" << endl;
        return false;
    }
    
    // Sort data by x-values in ascending order
    sortData();
    
    return true;
}

/**
 * Print first few data points for verification
 * @param numPoints - number of points to print
 */
void DataReader::print(size_t numPoints) const {
    cout << "\nFirst " << min(numPoints, xData.size()) << " data points:" << endl;
    for (size_t i = 0; i < min(numPoints, xData.size()); i++) {
        cout << "  x[" << i << "] = " << xData[i] 
                  << ", y[" << i << "] = " << yData[i] << endl;
    }
    cout << endl;
}

/**
 * Sort data by x-values in ascending order
 * For cubic spline interpolation
 * Keeps x,y pairs together during sorting
 */
void DataReader::sortData() {
    if (xData.size() != yData.size() || xData.empty()) {
        return;
    }
    
    // Check if already sorted
    if (isSorted()) {
        cout << "Data is already sorted in ascending order" << endl;
        return;
    }
    
    // Create vector of indices
    vector<size_t> indices(xData.size());
    for (size_t i = 0; i < indices.size(); i++) {
        indices[i] = i;
    }
    
    // Sort indices based on x-values
    sort(indices.begin(), indices.end(),
              [this](size_t i1, size_t i2) { return xData[i1] < xData[i2]; });
    
    // Create temporary copies
    vector<double> sortedX(xData.size());
    vector<double> sortedY(yData.size());
    
    // Reorder data using sorted indices
    for (size_t i = 0; i < indices.size(); i++) {
        sortedX[i] = xData[indices[i]];
        sortedY[i] = yData[indices[i]];
    }
    
    // Replace original data with sorted data
    xData = sortedX;
    yData = sortedY;
    
    cout << "Data sorted in ascending order by x-values" << endl;
}

/**
 * Check if data is sorted in ascending order
 * @return true if x-values are in ascending order
 */
bool DataReader::isSorted() const {
    for (size_t i = 1; i < xData.size(); i++) {
        if (xData[i] < xData[i-1]) {
            return false;
        }
    }
    return true;
}

/**
 * Find TMS peak (most positive x-value peak above baseline) and shift all data
 * so that the TMS peak maximum is at x = 0.0
 * 
 * @param baselineAdjustment - threshold for considering a point "above baseline"
 * @return the shift amount applied
 */
double DataReader::findAndShiftTMS(double baselineAdjustment) {
    if (xData.empty()) {
        cerr << "Error: No data to process for TMS shift" << endl;
        return 0.0;
    }
    
    // Find the most positive x-value that has y-value above baseline
    double tmsLocation = xData[0];
    double maxY = yData[0];
    
    // Look for peaks above baseline, tracking the rightmost (most positive x) one
    for (size_t i = 0; i < xData.size(); i++) {
        if (yData[i] > baselineAdjustment) {
            // Check if this is a local maximum
            bool isLocalMax = true;
            
            // Check neighbors
            if (i > 0 && yData[i] < yData[i-1]) {
                isLocalMax = false;
            }
            if (i < xData.size() - 1 && yData[i] < yData[i+1]) {
                isLocalMax = false;
            }
            
            // If it's a local max and more positive than current TMS, update
            if (isLocalMax && (xData[i] > tmsLocation || yData[i] > maxY)) {
                if (xData[i] >= tmsLocation) {  // Prioritize most positive x
                    tmsLocation = xData[i];
                    maxY = yData[i];
                }
            }
        }
    }
    
    // Shift all x-values so TMS is at 0.0
    double shift = tmsLocation;
    for (size_t i = 0; i < xData.size(); i++) {
        xData[i] -= shift;
    }
    
    cout << "TMS peak found at original x = " << tmsLocation << endl;
    cout << "Applied shift of " << shift << " ppm for TMS calibration" << endl;
    cout << endl;
    
    return shift;
}

// Legacy baseline correction, not implemented
double DataReader::correctBaseline() {
    // No baseline correction applied 
    return 0.0;
}