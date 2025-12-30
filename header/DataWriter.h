#ifndef DATAWRITER_H
#define DATAWRITER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

/**
 * DataWriter class - Utility for writing data to files at each processing step
 * 
 * This allows us to:
 * - Preserve original data
 * - Debug each transformation step
 * - Visualize intermediate results
 * - Verify correctness at each stage
 */
class DataWriter {
public:
    /**
     * Write x,y data pairs to a file
     * @param filename - output filename
     * @param xData - x values
     * @param yData - y values
     * @param header - optional header comment
     * @return true if successful
     */
    static bool writeData(const string& filename,
                         const vector<double>& xData,
                         const vector<double>& yData,
                         const string& header = "");
    
    /**
     * Write spline-evaluated data at many points for plotting
     * @param filename - output filename
     * @param spline - cubic spline to evaluate
     * @param xMin - minimum x value
     * @param xMax - maximum x value
     * @param numPoints - number of points to evaluate
     * @return true if successful
     */
    static bool writeSplineData(const string& filename,
                               const class CubicSpline& spline,
                               double xMin, double xMax,
                               int numPoints = 1000);
    
    /**
     * Write peak information to file for easy plotting
     * @param filename - output filename
     * @param peaks - vector of detected peaks
     * @param baseline - baseline level
     * @return true if successful
     */
    static bool writePeakData(const string& filename,
                             const vector<struct Peak>& peaks,
                             double baseline);
};

#endif // DATAWRITER_H
