#include "CubicSpline.h"
#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

/**
 * Constructor
 */
CubicSpline::CubicSpline() : computed(false) {
}

/**
 * Compute spline coefficients using natural cubic spline
 * 
 * Natural boundary conditions: second derivative = 0 at endpoints
 * Solves tridiagonal system for second derivatives M_i at each knot
 */
bool CubicSpline::compute(const vector<double>& xData, const vector<double>& yData) {
    if (xData.size() != yData.size() || xData.size() < 2) {
        cerr << "Error: Invalid data for spline computation" << endl;
        return false;
    }
    
    x = xData;
    y = yData;
    
    size_t n = x.size();
    b.resize(n);
    c.resize(n);
    d.resize(n);
    
    cout << "Computing natural cubic spline for " << n << " data points..." << endl;
    
    // Special case: only 2 points (linear interpolation)
    if (n == 2) {
        b[0] = (y[1] - y[0]) / (x[1] - x[0]);
        c[0] = 0.0;
        d[0] = 0.0;
        b[1] = b[0];
        c[1] = 0.0;
        d[1] = 0.0;
        computed = true;
        cout << "  Linear interpolation (2 points)" << endl;
        return true;
    }
    
    // Compute interval widths h_i = x_{i+1} - x_i
    vector<double> h(n-1);
    for (size_t i = 0; i < n-1; i++) {
        h[i] = x[i+1] - x[i];
        if (h[i] <= 0) {
            cerr << "Error: x values must be strictly increasing" << endl;
            return false;
        }
    }
    
    // Set up tridiagonal system for second derivatives M_i
    // Natural boundary conditions: M_0 = 0, M_{n-1} = 0
    // So we only solve for M_1, M_2, ..., M_{n-2}
    
    size_t m = n - 2;  // Number of unknowns
    
    if (m == 0) {
        // Only 2 points, already handled above
        computed = true;
        return true;
    }
    
    // Build tridiagonal matrix A and right-hand side vector rhs
    mat A = zeros<mat>(m, m);
    vec rhs = zeros<vec>(m);
    
    // First row (i=1)
    A(0, 0) = 2.0 * (h[0] + h[1]);
    if (m > 1) {
        A(0, 1) = h[1];
    }
    rhs(0) = 6.0 * ((y[2] - y[1]) / h[1] - (y[1] - y[0]) / h[0]);
    
    // Interior rows (i=2 to n-3)
    for (size_t i = 1; i < m-1; i++) {
        A(i, i-1) = h[i];
        A(i, i) = 2.0 * (h[i] + h[i+1]);
        A(i, i+1) = h[i+1];
        rhs(i) = 6.0 * ((y[i+2] - y[i+1]) / h[i+1] - (y[i+1] - y[i]) / h[i]);
    }
    
    // Last row (i=n-2)
    if (m > 1) {
        A(m-1, m-2) = h[m-1];
        A(m-1, m-1) = 2.0 * (h[m-1] + h[m]);
        rhs(m-1) = 6.0 * ((y[m+1] - y[m]) / h[m] - (y[m] - y[m-1]) / h[m-1]);
    }
    
    // Solve tridiagonal system using Armadillo
    vec M_interior = solve(A, rhs);
    
    // Construct full M vector with boundary conditions
    vector<double> M(n);
    M[0] = 0.0;  // Natural BC
    for (size_t i = 0; i < m; i++) {
        M[i+1] = M_interior(i);
    }
    M[n-1] = 0.0;  // Natural BC
    
    // Compute spline coefficients for each interval [x_i, x_{i+1}]
    // S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
    // where a_i = y_i (stored implicitly in evaluate)
    for (size_t i = 0; i < n-1; i++) {
        d[i] = (M[i+1] - M[i]) / (6.0 * h[i]);
        c[i] = M[i] / 2.0;
        b[i] = (y[i+1] - y[i]) / h[i] - h[i] * (2.0 * M[i] + M[i+1]) / 6.0;
    }
    
    // Set last point coefficients (not used in evaluation but kept for completeness)
    b[n-1] = b[n-2];
    c[n-1] = c[n-2];
    d[n-1] = d[n-2];
    
    cout << "  Tridiagonal system solved (" << m << " unknowns)" << endl;
    cout << "  Spline coefficients computed for " << (n-1) << " intervals" << endl;
    
    computed = true;
    return true;
}

/**
 * Evaluate spline at given x value
 * S_i(x) = y_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
 */
double CubicSpline::evaluate(double xVal) const {
    if (!computed || x.empty()) {
        return 0.0;
    }
    
    size_t n = x.size();
    
    // Find the interval [x_i, x_{i+1}] containing xVal
    size_t i = 0;
    
    // Handle extrapolation
    if (xVal <= x[0]) {
        i = 0;
    } else if (xVal >= x[n-1]) {
        i = n - 2;
    } else {
        // Binary search for efficiency (x is sorted)
        size_t left = 0, right = n - 1;
        while (right - left > 1) {
            size_t mid = (left + right) / 2;
            if (xVal < x[mid]) {
                right = mid;
            } else {
                left = mid;
            }
        }
        i = left;
    }
    
    // Evaluate spline polynomial for interval i
    double dx = xVal - x[i];
    return y[i] + b[i]*dx + c[i]*dx*dx + d[i]*dx*dx*dx;
}

/**
 * Evaluate spline derivative at given x value
 * S'_i(x) = b_i + 2*c_i*(x-x_i) + 3*d_i*(x-x_i)^2
 */
double CubicSpline::evaluateDerivative(double xVal) const {
    if (!computed || x.empty()) {
        return 0.0;
    }
    
    size_t n = x.size();
    
    // Find the interval [x_i, x_{i+1}] containing xVal
    size_t i = 0;
    
    if (xVal <= x[0]) {
        i = 0;
    } else if (xVal >= x[n-1]) {
        i = n - 2;
    } else {
        // Binary search
        size_t left = 0, right = n - 1;
        while (right - left > 1) {
            size_t mid = (left + right) / 2;
            if (xVal < x[mid]) {
                right = mid;
            } else {
                left = mid;
            }
        }
        i = left;
    }
    
    // Evaluate derivative of spline polynomial for interval i
    double dx = xVal - x[i];
    return b[i] + 2.0*c[i]*dx + 3.0*d[i]*dx*dx;
}

/**
 * Find crossings with horizontal line y = yVal
 * Uses sign change detection and bisection method refinement
 */
vector<double> CubicSpline::findCrossings(double yVal, double xMin, double xMax) const {
    vector<double> crossings;
    
    if (!computed || x.empty()) {
        return crossings;
    }
    
    // Sample the spline and look for sign changes
    int numSamples = 1000;
    double dx = (xMax - xMin) / numSamples;
    
    double prevVal = evaluate(xMin) - yVal;
    double prevX = xMin;
    
    for (int i = 1; i <= numSamples; i++) {
        double xSample = xMin + i * dx;
        double currVal = evaluate(xSample) - yVal;
        
        // Check for sign change
        if (prevVal * currVal < 0) {
            // Found a crossing, refine with bisection method
            double xLeft = prevX;
            double xRight = xSample;
            double fLeft = prevVal;
            double fRight = currVal;
            
            // Bisection iterations
            for (int iter = 0; iter < 50; iter++) {
                double xMid = (xLeft + xRight) / 2.0;
                double fMid = evaluate(xMid) - yVal;
                
                // Check convergence
                if (abs(fMid) < 1e-8 || abs(xRight - xLeft) < 1e-10) {
                    crossings.push_back(xMid);
                    break;
                }
                
                // Update interval based on sign
                if (fLeft * fMid < 0) {
                    // Root is in left half
                    xRight = xMid;
                    fRight = fMid;
                } else {
                    // Root is in right half
                    xLeft = xMid;
                    fLeft = fMid;
                }
                
                // If last iteration, add the midpoint
                if (iter == 49) {
                    crossings.push_back((xLeft + xRight) / 2.0);
                }
            }
        }
        
        prevVal = currVal;
        prevX = xSample;
    }
    
    cout << "  Found " << crossings.size() << " baseline crossings" << endl;
    
    return crossings;
}
