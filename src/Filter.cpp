#include "Filter.h"
#include <iostream>

using namespace std;

/**
 * Apply boxcar filter multiple times
 */
vector<double> Filter::applyBoxcar(const vector<double>& data, 
                                       int filterSize, int numPasses) {
    cout << "Applying " << numPasses << "-pass boxcar filter (size " 
              << filterSize << ")..." << endl;
    
    vector<double> result = data;
    for (int pass = 0; pass < numPasses; pass++) {
        result = boxcarPass(result, filterSize);
        cout << "  Pass " << (pass + 1) << " complete" << endl;
    }
    
    return result;
}

/**
 * Single pass of boxcar filter with cyclic boundary conditions
 * 
 * Formula: B(yi) = sum(y[i-k]...y[i]...y[i+k]) / n
 * where k = (n-1)/2 and boundaries wrap around cyclically
 */
vector<double> Filter::boxcarPass(const vector<double>& data, int filterSize) {
    if (filterSize <= 0 || data.empty()) {
        return data;
    }
    
    size_t n = data.size();
    vector<double> result(n);
    
    // Calculate half-width of filter window
    int halfWidth = (filterSize - 1) / 2;
    
    // Apply boxcar filter to each point
    for (size_t i = 0; i < n; i++) {
        double sum = 0.0;
        
        // Sum over filter window with reflection boundaries
        for (int j = -halfWidth; j <= halfWidth; j++) {
            // Calculate index with reflection at boundaries
            int idx = static_cast<int>(i) + j;
            
            // Handle negative indices (mirror at start)
            if (idx < 0) {
                idx = -idx;
            }
            // Handle indices past end (mirror at end)
            if (idx >= static_cast<int>(n)) {
                idx = 2 * static_cast<int>(n) - idx - 2;
            }
            // Safety bounds (for double reflection)
            if (idx < 0) {
                idx = 0;
            }
            if (idx >= static_cast<int>(n)) {
                idx = n - 1;
            }
            
            sum += data[idx];
        }
        
        // Calculate average
        result[i] = sum / filterSize;
    }
    
    return result;
}

/**
 * Apply Savitzky-Golay filter multiple times
 */
vector<double> Filter::applySavitzkyGolay(const vector<double>& data, 
                                              int filterSize, int numPasses) {
    cout << "Applying " << numPasses << "-pass Savitzky-Golay filter (size " 
              << filterSize << ")..." << endl;
    
    vector<double> result = data;
    for (int pass = 0; pass < numPasses; pass++) {
        result = sgPass(result, filterSize);
        cout << "  Pass " << (pass + 1) << " complete" << endl;
    }
    
    return result;
}

/**
 * Single pass of Savitzky-Golay filter
 * 
 * Uses quadratic polynomial fitting with pre-computed convolution coefficients
 * Coefficients from: Savitzky & Golay, Analytical Chemistry, 36, 1627 (1964)
 */
vector<double> Filter::sgPass(const vector<double>& data, int filterSize) {
    if (filterSize != 5 && filterSize != 11 && filterSize != 17) {
        cerr << "Warning: SG filter size should be 5, 11, or 17. Using 5." << endl;
        filterSize = 5;
    }
    
    if (data.empty()) {
        return data;
    }
    
    size_t n = data.size();
    vector<double> result(n);
    
    // Savitzky-Golay convolution coefficients for quadratic polynomial
    vector<double> coeffs;
    double norm;
    
    if (filterSize == 5) {
        // 5-point quadratic smoothing coefficients
        coeffs = {-3.0, 12.0, 17.0, 12.0, -3.0};
        norm = 35.0;
    } else if (filterSize == 11) {
        // 11-point quadratic smoothing coefficients
        coeffs = {-36.0, 9.0, 44.0, 69.0, 84.0, 89.0, 84.0, 69.0, 44.0, 9.0, -36.0};
        norm = 429.0;
    } else { // filterSize == 17
        // 17-point quadratic smoothing coefficients
        coeffs = {-21.0, -6.0, 7.0, 18.0, 27.0, 34.0, 39.0, 42.0, 43.0, 
                  42.0, 39.0, 34.0, 27.0, 18.0, 7.0, -6.0, -21.0};
        norm = 323.0;
    }
    
    int halfWidth = (filterSize - 1) / 2;
    
    // Apply SG filter to each point with reflection boundaries
    for (size_t i = 0; i < n; i++) {
        double sum = 0.0;
        
        // Convolve with coefficients
        for (int j = -halfWidth; j <= halfWidth; j++) {
            // Calculate index with reflection at boundaries
            int idx = static_cast<int>(i) + j;
            
            // Handle negative indices (mirror at start)
            if (idx < 0) {
                idx = -idx;
            }
            // Handle indices past end (mirror at end)
            if (idx >= static_cast<int>(n)) {
                idx = 2 * static_cast<int>(n) - idx - 2;
            }
            // Safety bounds (for double reflection)
            if (idx < 0) {
                idx = 0;
            }
            if (idx >= static_cast<int>(n)) {
                idx = n - 1;
            }
            
            // Apply convolution coefficient
            int coeffIdx = j + halfWidth;
            sum += coeffs[coeffIdx] * data[idx];
        }
        
        // Normalize
        result[i] = sum / norm;
    }
    
    return result;
}