#ifndef FILTER_H
#define FILTER_H

#include <vector>

using namespace std;

/**
 * Filter class - Applies smoothing filters to data
 * 
 * Supports:
 * - Boxcar (moving average) filter with cyclic boundary conditions
 * - Savitzky-Golay filter (5, 11, or 17 point)
 */
class Filter {
public:
    /**
     * Apply boxcar (moving average) filter
     * @param data - input data to filter
     * @param filterSize - size of filter window (must be odd)
     * @param numPasses - number of times to apply the filter
     * @return filtered data
     */
    static vector<double> applyBoxcar(const vector<double>& data, 
                                          int filterSize, int numPasses);
    
    /**
     * Apply Savitzky-Golay filter
     * @param data - input data to filter
     * @param filterSize - size of filter window (5, 11, or 17)
     * @param numPasses - number of times to apply the filter
     * @return filtered data
     */
    static vector<double> applySavitzkyGolay(const vector<double>& data, 
                                                  int filterSize, int numPasses);

private:
    // Helper function for single pass of boxcar filter
    static vector<double> boxcarPass(const vector<double>& data, int filterSize);
    
    // Helper function for single pass of SG filter
    static vector<double> sgPass(const vector<double>& data, int filterSize);
};

#endif // FILTER_H
