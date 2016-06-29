/*
 * links.h
 *
 *  Created on: 15.03.2010
 *      Author: daniel
 *
 * Classes for GLM link functions.
 *
 */

#ifndef LINKS_H_
#define LINKS_H_

#include <rcppExport.h>

static const double THRESH = 30.;
static const double MTHRESH = -30.;
static const double INVEPS = 1 / DOUBLE_EPS;


// ***************************************************************************************************//

// The general virtual Link class, which includes
// the link function itself,
// the inverse link function (or response function)
// and the derivative of the response function.
class Link
{
public:
    // the link function itself
    virtual double
    linkfun(double mu) const = 0;

    // the inverse link function (or response function)
    virtual double
    linkinv(double eta) const = 0;

    // and the derivative of the response function
    virtual double
    mu_eta(double eta) const = 0;

    // we need a virtual destructor here,
    // cf. Accelerated C++ pp. 242 ff.
    virtual ~Link(){}
};

// ***************************************************************************************************//

// the Logit link class
// inspiration taken from R/src/library/stats/src/family.c,
// because for numerical stability of the IWLS we REALLY NEED the thresholds...
class LogitLink : public Link
{
private:

    /**
     * Evaluate x/(1 - x). An inline function is used so that x is
     * evaluated once only.
     */
    static double x_d_omx(double x)
    {
        if (x < 0 || x > 1)
            Rf_error("Value %d out of range (0, 1)", x);
        return x / (1 - x);
    }

    /**
     * Evaluate x/(1 + x). An inline function is used so that x is
     * evaluated once only.
     */
    static double x_d_opx(double x)
    {
        return x / (1 + x);
    }

public:
    // logit
    double
    linkfun(double mu) const
    {
        return log(x_d_omx(mu));
    }

    // inverse logit
    double
    linkinv(double eta) const
    {
        double tmp = (eta < MTHRESH) ? DOUBLE_EPS :
                     ((eta > THRESH) ? INVEPS : exp(eta));
        return x_d_opx(tmp);
    }

    // function(eta) h'(eta) where h is the inverse logit
    double
    mu_eta(double eta) const
    {
        double opexp = 1 + exp(eta);
        double ret = (eta > THRESH || eta < MTHRESH) ? DOUBLE_EPS :
                      exp(eta) / (opexp * opexp);
        return ret;
    }
};


// ***************************************************************************************************//

// the Cauchit link class
class CauchitLink : public Link
{
public:
    // ctr
    CauchitLink() :
        thresh(- Rf_qcauchy(DOUBLE_EPS, 0.0, 1.0, 1, 0))
        {
        }

    // cauchit == cauchy quantile
    double
    linkfun(double mu) const
    {
        return Rf_qcauchy(mu, 0.0, 1.0, 1, 0);
    }

    // inverse cauchit == cauchy cdf
    // compare in R: binomial("cauchit")$linkinv
    double
    linkinv(double eta) const
    {
        eta = fmin(fmax(eta, - thresh), thresh);
        return Rf_pcauchy(eta, 0.0, 1.0, 1, 0);
    }

    // cauchy density
    double
    mu_eta(double eta) const
    {
        return fmax(Rf_dcauchy(eta, 0.0, 1.0, 0), DOUBLE_EPS);
    }

private:
    const double thresh;
};

// ***************************************************************************************************//

// the Probit link class
class ProbitLink : public Link
{
public:
    // ctr
    ProbitLink() :
        thresh(- Rf_qnorm5(DOUBLE_EPS, 0.0, 1.0, 1, 0))
        {
        }

    // probit == normal quantile
    double
    linkfun(double mu) const
    {
        return Rf_qnorm5(mu, 0.0, 1.0, 1, 0);
    }

    // inverse probit == normal cdf
    // compare in R: binomial("probit")$linkinv
    double
    linkinv(double eta) const
    {
        eta = fmin(fmax(eta, - thresh), thresh);
        return Rf_pnorm5(eta, 0.0, 1.0, 1, 0);
    }

    // normal density
    double
    mu_eta(double eta) const
    {
        return fmax(Rf_dnorm4(eta, 0.0, 1.0, 0), DOUBLE_EPS);
    }

private:
    const double thresh;
};

// ***************************************************************************************************//

// the complementary log-log link class
class CloglogLink : public Link
{
public:
    // cloglog
    double
    linkfun(double mu) const
    {
        return log(- log(1.0 - mu));
    }

    // compare in R: binomial("cloglog")$linkinv
    double
    linkinv(double eta) const
    {
        return fmax(fmin(- expm1(- exp(eta)), 1 - DOUBLE_EPS), DOUBLE_EPS);
    }

    // compare in R: binomial("cloglog")$mu.eta
    double
    mu_eta(double eta) const
    {
        eta = fmin(eta, 700.0);
        return fmax(exp(eta) * exp(- exp(eta)), DOUBLE_EPS);
    }
};

// ***************************************************************************************************//

// the inverse link class
class InverseLink : public Link
{
public:
    // inverse is the link
    double
    linkfun(double mu) const
    {
        return 1.0 / mu;
    }

    // and also the response function
    double
    linkinv(double eta) const
    {
        return 1.0 / eta;
    }

    // finally the response derivative
    double
    mu_eta(double eta) const
    {
       return - 1.0 / eta / eta;
    }
};

// ***************************************************************************************************//

// the log link class
class LogLink : public Link
{
public:
    // log is the link
    double
    linkfun(double mu) const
    {
        return log(mu);
    }

    // so the response function is exp
    double
    linkinv(double eta) const
    {
        return fmax(exp(eta), DOUBLE_EPS);
    }

    // finally the response derivative is also exp
    double
    mu_eta(double eta) const
    {
        return fmax(exp(eta), DOUBLE_EPS);
    }
};

// ***************************************************************************************************//

// the identity link class
class IdentityLink : public Link
{
public:
    // identity is the link ...
    double
    linkfun(double mu) const
    {
        return mu;
    }

    // ... and response function
    double
    linkinv(double eta) const
    {
        return eta;
    }

    // finally the response derivative is 1
    double
    mu_eta(double eta) const
    {
        return 1.0;
    }
};
#endif /* LINKS_H_ */
