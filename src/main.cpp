#include <iostream>
#include <chrono>
#include <fstream>
#include <algorithm>
#include "Config.h"
#include "DataReader.h"
#include "Filter.h"
#include "CubicSpline.h"
#include "PeakDetector.h"
#include "DataWriter.h"

using namespace std;
using namespace chrono;

/**
 * Main program for NMR Spectrum Analysis
 * 
 * CSC/MAT 335 - Project 1
 * Cubic Splines, Simple Filters, and Numerical Integration
 * 
 * This program:
 * 1. Reads configuration from nmr.in
 * 2. Reads NMR spectrum data
 * 3. Applies TMS calibration (shifts spectrum)
 * 4. Applies smoothing filter (boxcar or Savitzky-Golay)
 * 5. Fits natural cubic spline to filtered data
 * 6. Detects peaks above baseline
 * 7. Integrates peak areas using specified method
 * 8. Calculates relative hydrogen counts
 * 9. Outputs results to file
 */
int main(int argc, char* argv[]) {
    auto startTime = high_resolution_clock::now();
    
    // Determine config file name (default: nmr.in)
    string configFile = "nmr.in";
    if (argc > 1) {
        configFile = argv[1];
    }
    
    cout << "================================================" << endl;
    cout << "     NMR Spectrum Analysis Program" << endl;
    cout << "     CSC/MAT 335 - Project 1" << endl;
    cout << "================================================" << endl;
    cout << endl;
    
    // Read configuration
    Config config;
    if (!config.readFromFile(configFile)) {
        cerr << "Failed to read configuration file: " << configFile << endl;
        return 1;
    }
    config.print();
    
    // Read NMR data
    DataReader data;
    if (!data.readFromFile(config.inputFilename)) {
        cerr << "Failed to read data file: " << config.inputFilename << endl;
        return 1;
    }
    
    // Find TMS peak and shift spectrum
    double tmsShift = data.findAndShiftTMS(config.baselineAdjustment);
    
    // Legacy baseline correction, not implemented
    double baselineValue = data.correctBaseline();
    
    // Save shifted and baseline-corrected data
    DataWriter::writeData("shifted_data.txt", data.xData, data.yData, 
                         "Data after TMS calibration (shifted " + to_string(tmsShift) + 
                         " ppm) and baseline correction (baseline=" + to_string(baselineValue) + ")");

    cout << endl;
    
    // Apply filter (if enabled)
    vector<double> filteredY = data.yData;
    if (config.filterType == 1) {
        // Boxcar filter
        filteredY = Filter::applyBoxcar(data.yData, config.filterSize, config.filterPasses);
    } else if (config.filterType == 2) {
        // Savitzky-Golay filter
        filteredY = Filter::applySavitzkyGolay(data.yData, config.filterSize, config.filterPasses);
    } else {
        cout << "Filtering disabled (filter type = 0)" << endl;
    }
    
    // Save filtered data (if filtering was applied)
    if (config.filterType != 0) {
        DataWriter::writeData("filtered_data.txt", data.xData, filteredY,
                             "Data after " + config.getFilterTypeName() + " filtering");
    }
    cout << endl;
    
    // Fit cubic spline to (filtered) data
    CubicSpline spline;
    if (!spline.compute(data.xData, filteredY)) {
        cerr << "Failed to compute cubic spline" << endl;
        return 1;
    }
    
    // Save spline-evaluated data for visualization 
    if (spline.isComputed() && !data.xData.empty()) {
        double xMin = *std::min_element(data.xData.begin(), data.xData.end());
        double xMax = *std::max_element(data.xData.begin(), data.xData.end());
        DataWriter::writeSplineData("spline_fit.txt", spline, xMin, xMax, 2000);
    }
    cout << endl;
    
    // Detect peaks
    vector<Peak> peaks = PeakDetector::detectPeaks(spline, data.xData, filteredY, config.baselineAdjustment);
    cout << endl;
    
    // Integrate peaks
    PeakDetector::integratePeaks(peaks, spline, config.integrationType, config.tolerance);
    cout << endl;
    
    // Calculate hydrogen ratios
    PeakDetector::calculateHydrogens(peaks);
    
    // Save peak data for plotting/analysis
    DataWriter::writePeakData("peak_data.txt", peaks, config.baselineAdjustment);
    cout << endl;
    
    // Display results
    cout << "Techniques" << endl;
    cout << "===============================" << endl;
    cout << config.getIntegrationTypeName() << " Integration" << endl;
    cout << endl;
    
    cout << "Plot File Data" << endl;
    cout << "===============================" << endl;
    cout << "File: " << config.inputFilename << endl;
    cout << "Plot shifted " << tmsShift << " ppm for TMS calibration" << endl;
    cout << "Baseline corrected (subtracted " << baselineValue << ")" << endl;
    
    // Print formatted peak table to console
    PeakDetector::printPeaks(peaks);
    
    // Calculate execution time
    auto endTime = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(endTime - startTime);
    
    cout << "Analysis took " << (duration.count() / 1000.0) << " seconds." << endl;
    cout << endl;
    
    // Write results to output file
    ofstream outFile(config.outputFilename);
    if (outFile.is_open()) {
        outFile << "-=> NMR ANALYSIS <=-\n" << endl;
        outFile << "Program Options" << endl;
        outFile << "===============================" << endl;
        outFile << "Baseline Adjustment : " << config.baselineAdjustment << endl;
        outFile << "Tolerance           : " << config.tolerance << endl;
        outFile << "Filter Type         : " << config.getFilterTypeName() << endl;
        if (config.filterType != 0) {
            outFile << "Filter Size         : " << config.filterSize << endl;
            outFile << "Filter Passes       : " << config.filterPasses << endl;
        }
        outFile << "Integration Method  : " << config.getIntegrationTypeName() << endl;
        outFile << "\nTechniques" << endl;
        outFile << "===============================" << endl;
        outFile << config.getIntegrationTypeName() << " Integration" << endl;
        outFile << "\nPlot File Data" << endl;
        outFile << "===============================" << endl;
        outFile << "File: " << config.inputFilename << endl;
        outFile << "Plot shifted " << tmsShift << " ppm for TMS calibration" << endl;
        outFile << "Baseline corrected (subtracted " << baselineValue << ")\n" << endl;
        
        // Write peak table to file
        PeakDetector::printPeaks(outFile, peaks);
        
        outFile << "\nAnalysis took " << (duration.count() / 1000.0) << " seconds." << endl;
        
        outFile.close();
        cout << "Results written to: " << config.outputFilename << endl;
    } else {
        cerr << "Warning: Could not open output file: " << config.outputFilename << endl;
    }
    
    cout << "\n================================================" << endl;
    cout << "     Analysis Complete!" << endl;
    cout << "================================================" << endl;
    
    return 0;
}