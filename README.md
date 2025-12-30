# NMR Spectrum Analysis - Project 1
CSC 335 - Cubic Splines, Simple Filters, and Numerical Integration

## Project Overview
This program analyzes Nuclear Magnetic Resonance (NMR) spectroscopy data by:
1. Reading and calibrating spectral data (TMS peak alignment)
2. Applying smoothing filters to reduce noise
3. Fitting a natural cubic spline to the data
4. Detecting and locating peaks
5. Integrating peak areas using various numerical methods
6. Calculating relative hydrogen counts

## File Structure

### Source Files
- **main.cpp** - Main program that orchestrates the analysis workflow
- **Config.h/cpp** - Configuration file reader (reads nmr.in)
- **DataReader.h/cpp** - NMR data file reader and TMS calibration
- **Filter.h/cpp** - Data smoothing filters (boxcar and Savitzky-Golay)
- **CubicSpline.h/cpp** - Natural cubic spline fitting
- **Integration.h/cpp** - Numerical integration methods
- **PeakDetector.h/cpp** - Peak detection and analysis

### Build Files
- **Makefile** - Build automation

### Input Files
- **nmr.in** - Configuration file (see format below) //given test file
- **testdata.dat** - NMR spectrum data (x, y pairs) //given test file

### Output Files
- **analysis.txt** - Analysis results (configurable in nmr.in)

## Configuration File Format (nmr.in)

```
testdata.dat # Name of NMR input file
1000         # Baseline adjustment
1e-8         # Tolerance for numerical algorithms
1            # Type of Filter (0=none, 1=boxcar, 2=SG)
9            # Size of boxcar or SG filter (should be odd)
3            # Number of passes for the filter (ignored if Filter=0)
0            # Integration Technique (0=Newton-Cotes, 1=Romberg, 2=Adaptive, 3=Quadrature)
analysis.txt # Name of output file
```

### Parameters:
- **Line 1**: Input data filename
- **Line 2**: Baseline threshold for peak detection
- **Line 3**: Numerical tolerance for algorithms
- **Line 4**: Filter type (0=none, 1=boxcar, 2=Savitzky-Golay)
- **Line 5**: Filter window size (must be odd; 5, 11, or 17 for SG)
- **Line 6**: Number of filter passes
- **Line 7**: Integration method (0=Newton-Cotes, 1=Romberg, 2=Adaptive, 3=Gauss-Legendre)
- **Line 8**: Output filename

## Building and Running
The program may need to be ran from the data directory

### Compile the program:
```bash
make
```

### Run with default configuration (nmr.in):
```bash
./nmr_analysis
```

### Run with custom configuration file:
```bash
./nmr_analysis my_config.in
```

### Clean build artifacts:
```bash
make clean
```


## Algorithm Notes

### Natural Cubic Spline
- Requires solving a tridiagonal system of equations
- Natural boundary conditions: second derivative = 0 at endpoints

### Boxcar Filter
- Cyclic boundary conditions: reflect at edges
- More aggressive smoothing than Savitzky

### Savitzky-Golay Filter
- Polynomial smoothing filter
- Preserves peak shapes better than boxcar
- Coefficients from Savitzky & Golay (1964)

### Peak detection
- Uses midpoint to find peaks
- Chooses peaks from actual data, not spline

### Integration Methods
- All methods integrate the cubic spline (not raw data)
- Gauss-Legendre requires precomputed 64-point abscissas and weights

## Expected Output Format

```
-=> NMR ANALYSIS <=-

Program Options
===============================
Baseline Adjustment : 1000.0
Tolerance           : 1.0E-8
Filter Type         : Boxcar (Cyclic)
Filter Size         : 9
Filter Passes       : 3
Integration Method  : Newton-Cotes

Peak Begin           End              Location         Top              Area             Hydrogens
======= ================ ================ ================ ================ ================ =========
1       -4.880881966926  -4.449594398942  -4.567339208257  28098.892848     4.0135581393E+03 177
...

Analysis took 0.0025 seconds.
```
## Author
Cameron Lombardo
CSC 335 - Fall 2025