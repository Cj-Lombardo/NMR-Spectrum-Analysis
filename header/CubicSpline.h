#ifndef CUBICSPLINE_H
#define CUBICSPLINE_H

#include <vector>

using namespace std;

/**
 * CubicSpline class - Fits natural cubic spline to data
 * 
 * A natural cubic spline ensures:
 * - Spline passes through all data points
 * - First and second derivatives are continuous
 * - Second derivative is zero at endpoints (natural boundary conditions)
 * 
 * Requires solving a tridiagonal system of equations
 */
class CubicSpline {
private:
    vector<double> x;  // x data points
    vector<double> y;  // y data points
    vector<double> b;  // linear coefficients
    vector<double> c;  // quadratic coefficients
    vector<double> d;  // cubic coefficients
    bool computed;
    
public:
    CubicSpline();
    
    /**
     * Compute spline coefficients for given data
     * @param xData - x values (must be sorted)
     * @param yData - y values
     * @return true if successful
     */
    bool compute(const vector<double>& xData, const vector<double>& yData);
    
    /**
     * Evaluate spline at given x value
     * @param xVal - x value to evaluate at
     * @return interpolated y value
     */
    double evaluate(double xVal) const;
    
    /**
     * Evaluate spline derivative at given x value
     * @param xVal - x value to evaluate at
     * @return derivative value
     */
    double evaluateDerivative(double xVal) const;
    
    /**
     * Find x values where spline crosses a given y value
     * @param yVal - y value to find crossings for
     * @param xMin - minimum x to search
     * @param xMax - maximum x to search
     * @return vector of x values where spline crosses yVal
     */
    vector<double> findCrossings(double yVal, double xMin, double xMax) const;
    
    bool isComputed() const { return computed; }
};

#endif // CUBICSPLINE_H
