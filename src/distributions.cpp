/*
 * distributions.cpp
 *
 *  Created on: 15.03.2010
 *      Author: daniel
 *
 * Implementation of classes for exponential family distributions,
 * i.e. functions which should not be inlined by the compiler.
 *
 */

#include <distributions.h>
#include <string>
#include <types.h>
#include <rcppExport.h>

// ***************************************************************************************************//

// helper function for the binomial and Poisson loglikelihood
// see R/src/library/stats/src/family.c for the original R source code.
static inline double
y_log_y(double y, double mu)
{
    return (y > 0.0) ? (y * log(y / mu)) : 0.0;
}

// ***************************************************************************************************//


// binomial loglikelihood
// compare in R: binomial()$dev.resids
double
Binomial::loglik(const double *means) const
{
    double ret = 0.0;

    for(PosInt i = 0; i < responses.n_elem; ++i, ++means)
    {
        ret += weights(i) * (y_log_y(responses(i), *means) + y_log_y(1.0 - responses(i), 1.0 - *means));
    }

    return - ret;
}

// ***************************************************************************************************//


// Gaussian loglikelihood
// compare in R: poisson()$dev.resids
double
Gaussian::loglik(const double *means) const
{
    double ret = 0.0;

    for(PosInt i = 0; i < responses.n_elem; ++i, ++means)
    {
        ret += weights(i) * (responses(i) - *means) * (responses(i) - *means);
    }

    return - 0.5 * ret / phi;
}

// ***************************************************************************************************//


// Poisson loglikelihood
// compare in R: poisson()$dev.resids
double
Poisson::loglik(const double *means) const
{
    double ret = 0.0;

    for(PosInt i = 0; i < responses.n_elem; ++i, ++means)
    {
        ret += weights(i) * (responses(i) - *means - y_log_y(responses(i), *means));
    }

    return ret;
}
