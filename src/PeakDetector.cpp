#include "PeakDetector.h"
#include "Integration.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;

/**
 * Detect peaks in the spectrum
 * 
 * Algorithm:
 * 1. Find all baseline crossings using spline and bisection
 * 2. Pair crossings into peak regions (up-crossing to down-crossing)
 * 3. Peak location = midpoint (xa + xb)/2 as specified
 * 4. Find maximum from actual filtered data points (not spline)
 */
vector<Peak> PeakDetector::detectPeaks(const CubicSpline& spline,
                                            const vector<double>& xData,
                                            const vector<double>& yData,
                                            double baseline) {
    vector<Peak> peaks;
    
    if (!spline.isComputed() || xData.empty() || yData.empty()) {
        cerr << "Error: Invalid data for peak detection" << endl;
        return peaks;
    }
    
    if (xData.size() != yData.size()) {
        cerr << "Error: xData and yData size mismatch" << endl;
        return peaks;
    }
    
    cout << "Detecting peaks above baseline " << baseline << "..." << endl;
    
    // Get x-range from data
    double xMin = *min_element(xData.begin(), xData.end());
    double xMax = *max_element(xData.begin(), xData.end());
    
    // Find all baseline crossings using spline and bisection
    vector<double> crossings = spline.findCrossings(baseline, xMin, xMax);
    
    // Fix: Check for peaks the extend past data
    // Left edge: if spectrum starts above baseline
    double yAtMin = spline.evaluate(xMin);
    if (yAtMin > baseline && !crossings.empty()) {
        // Spectrum starts above baseline - add xMin as a virtual crossing
        crossings.insert(crossings.begin(), xMin);
    }
    
    // Check for peaks the extend past data
    // Right edge: if spectrum ends above baseline  
    double yAtMax = spline.evaluate(xMax);
    if (yAtMax > baseline && !crossings.empty()) {
        // Spectrum ends above baseline - add xMax as a virtual crossing
        crossings.push_back(xMax);
    }
    
    if (crossings.size() < 2) {
        cout << "  No complete peaks found (need at least 2 crossings)" << endl;
        return peaks;
    }
    
    // Sort crossings (should already be sorted, but ensure it)
    sort(crossings.begin(), crossings.end());
    
    cout << "  Found " << crossings.size() << " baseline crossings" << endl;
    
    // Group crossings into peak regions
    // A peak is between an up-crossing and the next down-crossing
    for (size_t i = 0; i < crossings.size() - 1; i++) {
        double xBegin = crossings[i];
        double xEnd = crossings[i + 1];
        
        // Check if there's actually a peak between these crossings
        // Sample a few points to verify it goes above baseline
        double xMid = (xBegin + xEnd) / 2.0;
        double yMid = spline.evaluate(xMid);
        
        if (yMid <= baseline) {
            // Not a peak, probably a valley - skip this pair
            continue;
        }
        
        // Find the maximum from ACTUAL DATA POINTS between xBegin and xEnd
        // (Not from spline - use the filtered data)
        double maxY = baseline;
        double maxX = xMid;  
        bool foundDataPoint = false;
        
        for (size_t j = 0; j < xData.size(); j++) {
            if (xData[j] >= xBegin && xData[j] <= xEnd) {
                if (yData[j] > maxY) {
                    maxY = yData[j];
                    maxX = xData[j];
                    foundDataPoint = true;
                }
            }
        }
        
        // If no data points found in range, skip this peak
        if (!foundDataPoint) {
            continue;
        }
        
        // Create peak structure
        Peak peak;
        peak.begin = xBegin;
        peak.end = xEnd;
        peak.location = (xBegin + xEnd) / 2.0;  // Midpoint as specified
        peak.maximum = maxY;  // From actual data point
        peak.area = 0.0;  // Will be calculated by integration
        peak.hydrogens = 0;  // Will be calculated later
        
        // Skip peak at 0 fix
        if (fabs(peak.location) < 0.02) {  // Within 0.02 ppm of zero
            continue;  // Skip this peak
        }
        
        peaks.push_back(peak);
    }
    
    cout << "  Found " << peaks.size() << " peaks" << endl;
    
    return peaks;
}

/**
 * Integrate all peaks using specified method
 */
void PeakDetector::integratePeaks(vector<Peak>& peaks,
                                 const CubicSpline& spline,
                                 int integrationType,
                                 double tolerance) {
    cout << "Integrating peaks..." << endl;
    
    for (size_t i = 0; i < peaks.size(); i++) {
        Peak& peak = peaks[i];
        
        // Integrate based on method
        switch (integrationType) {
            case 0:
                peak.area = Integration::newtonCotes(spline, peak.begin, peak.end, tolerance);
                break;
            case 1:
                peak.area = Integration::romberg(spline, peak.begin, peak.end, tolerance);
                break;
            case 2:
                peak.area = Integration::adaptive(spline, peak.begin, peak.end, tolerance);
                break;
            case 3:
                peak.area = Integration::gaussLegendre(spline, peak.begin, peak.end);
                break;
            default:
                cerr << "Unknown integration type: " << integrationType << endl;
                peak.area = 0.0;
        }
    }
}

/**
 * Calculate relative hydrogen counts
 */
void PeakDetector::calculateHydrogens(vector<Peak>& peaks) {
    if (peaks.empty()) {
        return;
    }
    
    // Find minimum area (represents 1 hydrogen)
    double minArea = peaks[0].area;
    for (const auto& peak : peaks) {
        if (peak.area > 0 && peak.area < minArea) {
            minArea = peak.area;
        }
    }
    
    // Calculate relative hydrogens
    for (auto& peak : peaks) {
        peak.hydrogens = static_cast<int>(round(peak.area / minArea));
    }
    
    cout << "Calculated hydrogen ratios (smallest peak = 1 H)" << endl;
}

/**
 * Print peaks in formatted table to output stream
 */
void PeakDetector::printPeaks(ostream& os, const vector<Peak>& peaks) {
    os << "\n";
    os << setw(7) << "Peak" << " "
              << setw(16) << "Begin" << " "
              << setw(16) << "End" << " "
              << setw(16) << "Location" << " "
              << setw(16) << "Top" << " "
              << setw(16) << "Area" << " "
              << setw(9) << "Hydrogens" << endl;
    
    os << string(7, '=') << " "
              << string(16, '=') << " "
              << string(16, '=') << " "
              << string(16, '=') << " "
              << string(16, '=') << " "
              << string(16, '=') << " "
              << string(9, '=') << endl;
    
    for (size_t i = 0; i < peaks.size(); i++) {
        const Peak& peak = peaks[i];
        os << setw(7) << (i+1) << " "
                  << setw(16) << fixed << setprecision(12) << peak.begin << " "
                  << setw(16) << fixed << setprecision(12) << peak.end << " "
                  << setw(16) << fixed << setprecision(12) << peak.location << " "
                  << setw(16) << fixed << setprecision(6) << peak.maximum << " "
                  << setw(16) << scientific << setprecision(10) << peak.area << " "
                  << setw(9) << peak.hydrogens << endl;
    }
    
    os << endl;
}

/**
 * Print peaks in formatted table to cout
 */
void PeakDetector::printPeaks(const vector<Peak>& peaks) {
    printPeaks(cout, peaks);
}