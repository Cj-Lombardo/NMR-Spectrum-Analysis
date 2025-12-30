#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "CubicSpline.h"

using namespace std;

/**
 * Integration class - Numerical integration methods
 * 
 * Supports:
 * - Newton-Cotes (composite trapezoidal/Simpson's)
 * - Romberg integration
 * - Adaptive quadrature
 * - Gauss-Legendre quadrature (64 points)
 */
class Integration {
public:
    /**
     * Integrate using composite Newton-Cotes (Simpson's rule)
     * @param spline - cubic spline to integrate
     * @param a - lower bound
     * @param b - upper bound
     * @param tolerance - error tolerance
     * @return integral value
     */
    static double newtonCotes(const CubicSpline& spline, double a, double b, double tolerance);
    
    /**
     * Integrate using Romberg method
     * @param spline - cubic spline to integrate
     * @param a - lower bound
     * @param b - upper bound
     * @param tolerance - error tolerance
     * @return integral value
     */
    static double romberg(const CubicSpline& spline, double a, double b, double tolerance);
    
    /**
     * Integrate using adaptive quadrature
     * @param spline - cubic spline to integrate
     * @param a - lower bound
     * @param b - upper bound
     * @param tolerance - error tolerance
     * @return integral value
     */
    static double adaptive(const CubicSpline& spline, double a, double b, double tolerance);
    
    /**
     * Integrate using Gauss-Legendre quadrature (64 points)
     * @param spline - cubic spline to integrate
     * @param a - lower bound
     * @param b - upper bound
     * @return integral value
     */
    static double gaussLegendre(const CubicSpline& spline, double a, double b);

private:
    // Helper function for trapezoidal rule
    static double trapezoid(const CubicSpline& spline, double a, double b, int n);
    
    // Helper for adaptive recursion
    static double adaptiveHelper(const CubicSpline& spline, double a, double b, 
                                 double tolerance, double fa, double fb, double fmid);
};

#endif // INTEGRATION_H
