/*
 * gpriors.cpp
 *
 *  Created on: 22.11.2012
 *      Author: daniel
 */


#include <gpriors.h>


// compute M(a, b)
static double
incInvGammaLogNormConst(double a, double b)
{
    double logNormConst;
    if(b > 0)
    {
        logNormConst = a * log(b) - Rf_pgamma(b, a, 1.0, 1, 1) - Rf_lgammafn(a);
    } else {
        logNormConst = log(a);
    }
    return logNormConst;
}

// Log prior density
double
IncInvGammaGPrior::logDens(double g) const
{
    return (- (a + 1.0) * log1p(g) - b / (g + 1.0) + incInvGammaLogNormConst(a, b));
}


// for this class we have a closed form for the log marginal likelihood
// resulting from the TBF approach
double
IncInvGammaGPrior::getTBFLogMargLik(double residualDeviance, int df) const
{
    return (incInvGammaLogNormConst(a, b) -
            incInvGammaLogNormConst(a + df / 2.0, b + residualDeviance / 2.0) +
            residualDeviance / 2.0);
}
