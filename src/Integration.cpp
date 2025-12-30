#include "Integration.h"
#include <iostream>
#include <cmath>

using namespace std;

/**
 * Newton-Cotes integration (composite Simpson's rule)
 * Uses adaptive subdivision until tolerance is met
 */
double Integration::newtonCotes(const CubicSpline& spline, double a, double b, double tolerance) {
    if (!spline.isComputed()) {
        cerr << "Error: Spline not computed" << endl;
        return 0.0;
    }
    
    // Start with a small number of intervals and double until convergence
    int n = 2;  // Start with 2 intervals (must be even for Simpson's)
    double prevIntegral = 0.0;
    double integral = 0.0;
    
    const int maxIterations = 20;
    
    for (int iter = 0; iter < maxIterations; iter++) {
        // Composite Simpson's rule: I = (h/3)[f(x0) + 4f(x1) + 2f(x2) + 4f(x3) + ... + f(xn)]
        double h = (b - a) / n;
        double sum = spline.evaluate(a) + spline.evaluate(b);
        
        // Odd indices get coefficient 4
        for (int i = 1; i < n; i += 2) {
            sum += 4.0 * spline.evaluate(a + i * h);
        }
        
        // Even indices (except endpoints) get coefficient 2
        for (int i = 2; i < n; i += 2) {
            sum += 2.0 * spline.evaluate(a + i * h);
        }
        
        integral = (h / 3.0) * sum;
        
        // Check convergence
        if (iter > 0 && abs(integral - prevIntegral) < tolerance) {
            return integral;
        }
        
        prevIntegral = integral;
        n *= 2;  // Double the number of intervals
    }
    
    // Return best estimate even if not fully converged
    return integral;
}

/**
 * Romberg integration
 * Uses Richardson extrapolation on trapezoidal rule
 */
double Integration::romberg(const CubicSpline& spline, double a, double b, double tolerance) {
    if (!spline.isComputed()) {
        cerr << "Error: Spline not computed" << endl;
        return 0.0;
    }
    
    const int maxLevel = 15;  // Maximum Romberg levels
    double R[maxLevel][maxLevel];
    
    // First column: trapezoidal rule with increasing subdivisions
    for (int i = 0; i < maxLevel; i++) {
        int n = static_cast<int>(pow(2, i));
        R[i][0] = trapezoid(spline, a, b, n);
        
        // Richardson extrapolation for higher order approximations
        for (int j = 1; j <= i; j++) {
            double factor = pow(4, j);
            R[i][j] = (factor * R[i][j-1] - R[i-1][j-1]) / (factor - 1.0);
        }
        
        // Check convergence (compare diagonal elements)
        if (i > 0 && abs(R[i][i] - R[i-1][i-1]) < tolerance) {
            return R[i][i];
        }
    }
    
    // Return best estimate
    return R[maxLevel-1][maxLevel-1];
}

/**
 * Adaptive quadrature
 * Uses recursive Simpson's rule with automatic subdivision
 */
double Integration::adaptive(const CubicSpline& spline, double a, double b, double tolerance) {
    if (!spline.isComputed()) {
        cerr << "Error: Spline not computed" << endl;
        return 0.0;
    }
    
    // Evaluate function at endpoints and midpoint
    double fa = spline.evaluate(a);
    double fb = spline.evaluate(b);
    double fmid = spline.evaluate((a + b) / 2.0);
    
    return adaptiveHelper(spline, a, b, tolerance, fa, fb, fmid);
}

/**
 * Gauss-Legendre quadrature (64 points)
 */
double Integration::gaussLegendre(const CubicSpline& spline, double a, double b) {
    if (!spline.isComputed()) {
        cerr << "Error: Spline not computed" << endl;
        return 0.0;
    }
    
    // 64-point Gauss-Legendre abscissas and weights
    // These are for the standard interval [-1, 1]
    // Only store positive values; use symmetry for negative
    const int n = 32;  // We'll use symmetry, so only need half
    
    // Nodes and weights for 64-point Gauss-Legendre quadrature (first 32)
    static const double nodes[32] = {
        0.0243502926634244325089558,
        0.0729931217877990394495429,
        0.1214628192961205544703765,
        0.1696444204239928180373136,
        0.2174236437400070841496487,
        0.2646871622087674163739642,
        0.3113228719902109561575127,
        0.3572201583376681159504426,
        0.4022701579639916036957668,
        0.4463660172534640879849477,
        0.4894031457070529574785263,
        0.5312794640198945456580139,
        0.5718956462026340342838781,
        0.6111553551723932502488530,
        0.6489654712546573398577612,
        0.6852363130542332425635584,
        0.7198818501716108268489402,
        0.7528199072605318966118638,
        0.7839723589433414076102205,
        0.8132653151227975597419233,
        0.8406292962525803627516915,
        0.8659993981540928197607834,
        0.8893154459951141058534040,
        0.9105221370785028057563807,
        0.9295691721319395758214902,
        0.9464113748584028160624815,
        0.9610087996520537189186141,
        0.9733268277899109637418535,
        0.9833362538846259569312993,
        0.9910133714767443207393824,
        0.9963401167719552793469245,
        0.9993050417357721394569056
    };
    
    static const double weights[32] = {
        0.0486909570091397203833654,
        0.0485754674415034269347991,
        0.0483447622348029571697695,
        0.0479993885964583077281262,
        0.0475401657148303086622822,
        0.0469681828162100173253263,
        0.0462847965813144172959532,
        0.0454916279274181444797710,
        0.0445905581637565630601347,
        0.0435837245293234533768279,
        0.0424735151236535890073398,
        0.0412625632426235286101563,
        0.0399537411327203413866569,
        0.0385501531786156291289625,
        0.0370551285402400460404151,
        0.0354722132568823838106931,
        0.0338051618371416093915655,
        0.0320579283548515535854675,
        0.0302346570724024788679741,
        0.0283396726142594832275113,
        0.0263774697150546586716918,
        0.0243527025687108733381776,
        0.0222701738083832541592983,
        0.0201348231535302093723403,
        0.0179517157756973430850453,
        0.0157260304760247193219660,
        0.0134630478967186425980608,
        0.0111681394601311288185905,
        0.0088467598263639477230309,
        0.0065044579689783628561174,
        0.0041470332605624676352875,
        0.0017832807216964329472961
    };
    
    // Transform from [-1, 1] to [a, b]
    double midpoint = (a + b) / 2.0;
    double halfwidth = (b - a) / 2.0;
    
    double sum = 0.0;
    
    // Use symmetry: integrate from -1 to 1 using both positive and negative nodes
    for (int i = 0; i < n; i++) {
        double x_pos = midpoint + halfwidth * nodes[i];
        double x_neg = midpoint - halfwidth * nodes[i];
        
        sum += weights[i] * (spline.evaluate(x_pos) + spline.evaluate(x_neg));
    }
    
    return halfwidth * sum;
}

/**
 * Helper: Trapezoidal rule
 */
double Integration::trapezoid(const CubicSpline& spline, double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (spline.evaluate(a) + spline.evaluate(b));
    
    for (int i = 1; i < n; i++) {
        sum += spline.evaluate(a + i * h);
    }
    
    return h * sum;
}

/**
 * Helper: Adaptive recursion using Simpson's rule
 * Recursively subdivides intervals until tolerance is met
 */
double Integration::adaptiveHelper(const CubicSpline& spline, double a, double b,
                                  double tolerance, double fa, double fb, double fmid) {
    double mid = (a + b) / 2.0;
    double h = b - a;
    
    // Simpson's rule on whole interval [a, b]
    double S_whole = (h / 6.0) * (fa + 4.0 * fmid + fb);
    
    // Simpson's rule on left half [a, mid]
    double leftMid = (a + mid) / 2.0;
    double f_leftMid = spline.evaluate(leftMid);
    double S_left = (h / 12.0) * (fa + 4.0 * f_leftMid + fmid);
    
    // Simpson's rule on right half [mid, b]
    double rightMid = (mid + b) / 2.0;
    double f_rightMid = spline.evaluate(rightMid);
    double S_right = (h / 12.0) * (fmid + 4.0 * f_rightMid + fb);
    
    double S_split = S_left + S_right;
    
    // Error estimate using Richardson extrapolation
    double error = abs(S_split - S_whole) / 15.0;
    
    // If error is acceptable, return refined estimate
    if (error < tolerance) {
        return S_split + (S_split - S_whole) / 15.0;  // Richardson correction
    }
    
    // Otherwise, recursively subdivide with tighter tolerance
    double left_integral = adaptiveHelper(spline, a, mid, tolerance / 2.0, 
                                         fa, fmid, f_leftMid);
    double right_integral = adaptiveHelper(spline, mid, b, tolerance / 2.0, 
                                          fmid, fb, f_rightMid);
    
    return left_integral + right_integral;
}
