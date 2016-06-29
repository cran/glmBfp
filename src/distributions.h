/*
 * distributions.h
 *
 *  Created on: 15.03.2010
 *      Author: daniel
 *
 * Classes for exponential family distributions.
 *
 */

#ifndef DISTRIBUTIONS_H_
#define DISTRIBUTIONS_H_

#include <types.h>
#include <rcppExport.h>


// ***************************************************************************************************//

// General virtual class for an exponential family distribution, including
// the variance function of the modelled mean mu
// and the loglikelihood of a whole vector of modelled means.
// We also always save the vector of scalar responses and weights, because naturally
// the loglikelihood of "means" depends on this!
class Distribution
{
public:
    // ctr:
    Distribution(const AVector& responses,
                 const AVector& weights) :
        responses(responses),
        weights(weights)
        {
        }

    // variance function
    virtual double
    variance(double mu) const = 0;

    // loglikelihood
    virtual double
    loglik(const double *means) const = 0;

    // we need a virtual destructor here,
    // cf. Accelerated C++ pp. 242 ff.
    virtual ~Distribution(){}

protected:
    // the vector of responses
    const AVector responses;

    // the vector of weights
    const AVector weights;
};

// ***************************************************************************************************//

// The binomial distribution
class Binomial : public Distribution
{
public:
    // ctr
    Binomial(const AVector& responses,
             const AVector& weights) :
        Distribution(responses,
                     weights)
        {
        }

    // variance function
    double
    variance(double mu) const
    {
        return mu * (1.0 - mu);
    }

    // loglikelihood
    double
    loglik(const double *means) const;
};

// ***************************************************************************************************//

// The Gaussian distribution
class Gaussian : public Distribution
{
public:
    // ctr
    Gaussian(const AVector& responses,
             const AVector& weights,
             double phi) :
        Distribution(responses, weights),
        phi(phi)
    {
    }

    // variance function
    double
    variance(double mu) const
    {
        return 1.0;
    }

    // loglikelihood
    double
    loglik(const double *means) const;

private:

    // the dispersion factor (variance sigma^2)
    const double phi;
};

// ***************************************************************************************************//

// The Poisson distribution
class Poisson : public Distribution
{
public:
    // ctr
    Poisson(const AVector& responses,
            const AVector& weights) :
        Distribution(responses,
                     weights)
    {
    }

    // variance function
    double
    variance(double mu) const
    {
        return mu;
    }

    // loglikelihood
    double
    loglik(const double *means) const;
};


#endif /* DISTRIBUTIONS_H_ */
