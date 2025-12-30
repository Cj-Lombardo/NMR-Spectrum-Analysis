#ifndef PEAKDETECTOR_H
#define PEAKDETECTOR_H

#include <ostream>
#include <vector>
#include "CubicSpline.h"


using namespace std;

/**
 * Peak structure - stores information about a detected peak
 */
struct Peak {
    double begin;      // x-value where peak starts (baseline crossing)
    double end;        // x-value where peak ends (baseline crossing)
    double location;   // x-value of peak maximum (midpoint between crossings)
    double maximum;    // y-value at peak maximum
    double area;       // integrated area of peak
    int hydrogens;     // relative number of hydrogens
};

/**
 * PeakDetector class - Finds and analyzes peaks in NMR spectrum
 */
class PeakDetector {
public:
    /**
     * Detect peaks in spline above baseline
     * @param spline - cubic spline fitted to data
     * @param xData - original x data points
     * @param yData - filtered y data points (for finding maximum)
     * @param baseline - baseline threshold
     * @return vector of detected peaks
     */
    static vector<Peak> detectPeaks(const CubicSpline& spline,
                                         const vector<double>& xData,
                                         const vector<double>& yData,
                                         double baseline);
    
    /**
     * Integrate peak areas using specified method
     * @param peaks - peaks to integrate
     * @param spline - cubic spline to integrate
     * @param integrationType - integration method (0-3)
     * @param tolerance - integration tolerance
     */
    static void integratePeaks(vector<Peak>& peaks,
                              const CubicSpline& spline,
                              int integrationType,
                              double tolerance);
    
    /**
     * Calculate relative hydrogen counts for peaks
     * @param peaks - peaks to analyze
     */
    static void calculateHydrogens(vector<Peak>& peaks);
    
    /**
     * Print peak analysis results (to cout)
     * @param peaks - peaks to print
     */
    static void printPeaks(const vector<Peak>& peaks);

    /**
     * Print peak analysis results to output stream
     * @param os - output stream
     * @param peaks - peaks to print
     */
    static void printPeaks(ostream& os, const vector<Peak>& peaks);
};

#endif // PEAKDETECTOR_H