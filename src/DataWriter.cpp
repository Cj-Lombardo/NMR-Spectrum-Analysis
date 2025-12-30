#include "DataWriter.h"
#include "CubicSpline.h"
#include "PeakDetector.h"
#include <iomanip>

using namespace std;

/**
 * Write x,y data to file
 */
bool DataWriter::writeData(const string& filename,
                          const vector<double>& xData,
                          const vector<double>& yData,
                          const string& header) {
    if (xData.size() != yData.size()) {
        cerr << "Error: x and y data sizes don't match" << endl;
        return false;
    }
    
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        cerr << "Error: Cannot open output file: " << filename << endl;
        return false;
    }
    
    // Write header if provided
    if (!header.empty()) {
        outFile << "# " << header << endl;
    }
    
    // Write data
    outFile << fixed << setprecision(6);
    for (size_t i = 0; i < xData.size(); i++) {
        outFile << xData[i] << " " << yData[i] << endl;
    }
    
    outFile.close();
    cout << "Data written to: " << filename << " (" << xData.size() << " points)" << endl;
    return true;
}

/**
 * Write spline-evaluated data
 */
bool DataWriter::writeSplineData(const string& filename,
                                const CubicSpline& spline,
                                double xMin, double xMax,
                                int numPoints) {
    if (!spline.isComputed()) {
        cerr << "Error: Spline not computed" << endl;
        return false;
    }
    
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        cerr << "Error: Cannot open output file: " << filename << endl;
        return false;
    }
    
    outFile << "# Cubic spline evaluated at " << numPoints << " points" << endl;
    outFile << fixed << setprecision(6);
    
    double dx = (xMax - xMin) / (numPoints - 1);
    for (int i = 0; i < numPoints; i++) {
        double x = xMin + i * dx;
        double y = spline.evaluate(x);
        outFile << x << " " << y << endl;
    }
    
    outFile.close();
    cout << "Spline data written to: " << filename << " (" << numPoints << " points)" << endl;
    return true;
}

/**
 * Write peak data for visualization
 */
bool DataWriter::writePeakData(const string& filename,
                              const vector<Peak>& peaks,
                              double baseline) {
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        cerr << "Error: Cannot open output file: " << filename << endl;
        return false;
    }
    
    outFile << "# Peak data for plotting" << endl;
    outFile << "# Format: peak_number, begin, end, location, maximum, area, hydrogens" << endl;
    outFile << "# Baseline: " << baseline << endl;
    outFile << fixed << setprecision(12);
    
    for (size_t i = 0; i < peaks.size(); i++) {
        outFile << (i+1) << " "
                << peaks[i].begin << " "
                << peaks[i].end << " "
                << peaks[i].location << " "
                << peaks[i].maximum << " "
                << scientific << peaks[i].area << " "
                << fixed << peaks[i].hydrogens << endl;
    }
    
    outFile.close();
    cout << "Peak data written to: " << filename << " (" << peaks.size() << " peaks)" << endl;
    return true;
}
