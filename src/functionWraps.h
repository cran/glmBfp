/*
 * functionWraps.h
 *
 *  Created on: 06.11.2009
 *      Author: daniel
 */

#ifndef FUNCTIONWRAPS_H_
#define FUNCTIONWRAPS_H_

#include <rcppExport.h>
#include <types.h>

// ***************************************************************************************************//

// wrap a simple R function,
// taking a double argument and returning a double.
class RFunction
{
public:
    // constructor
    RFunction(SEXP R_fun) :
        fun(R_fun)
        {}

//    // copy constructor
//    RFunction(const RFunction& old);
//
//    // destructor
//    ~RFunction();
//
//    // assignment operator
//    RFunction& operator=(const RFunction& rhs);

    // for use as a function (function object)
    double
    operator()(double x) const;

private:
    Rcpp::Function fun;
};

// ***************************************************************************************************//

// now wrap an R function taking vectors of fixed length
class VectorRFunction
{
public:
    // constructor
    VectorRFunction(SEXP R_function, R_len_t argsize) :
        fun(R_function),
        argsize(argsize)
        {}

    // for use as a function (function object)
    double
    operator()(const double* vector) const;

private:
    Rcpp::Function fun;
    PosInt argsize;
};

// ***************************************************************************************************//

// simple function value cache
class Cache
{
public:
    // save one pair of argument and function value
    void
    save(double arg, double val);

    // extract the arguments
    MyDoubleVector
    getArgs() const
    {
        return args;
    }

    // extract the function values
    MyDoubleVector
    getVals() const
    {
        return vals;
    }

    // clear the cache
    void
    clear();

    // look for an already computed function value
    double
    getValue(double arg) const;

    // initialize from an R list
    Cache(Rcpp::List& rcpp_list);

    // default ctr is allowed as well
    Cache(){};

    // convert to an R list
    Rcpp::List
    convert2list() const;

private:
    MyDoubleVector args;
    MyDoubleVector vals;
};

// ***************************************************************************************************//

// cache all function results
template<class Fun>
    class CachedFunction
    {
    public:
        // constructor takes function object reference and stores it.
        CachedFunction(Fun& function) :
            function(function)
        {
        }

        // return the cached values
        Cache
        getCache() const
        {
            return cache;
        }

        // for use as a function
        double
        operator()(double x);

    private:
        Fun& function;
        Cache cache;
    };

// ***************************************************************************************************//

// get derivative of a function
template<class Fun>
    class NumericDerivative
    {
    public:
        // save the function (object) reference
        NumericDerivative(Fun& function) :
            function(function)
        {
        }

        // for use as a function
        double
        operator()(double x) const;

    private:
        Fun& function;
    };

// ***************************************************************************************************//

// get more accurate derivative of a function
template<class Fun>
    class AccurateNumericDerivative
    {
    public:
        // save the function (object) reference
        AccurateNumericDerivative(Fun& function,
                                  double eps = 0.1) :
            function(function),
            h(eps),
            ntab(10),
            con(1.4),
            con2(con * con),
            big(DBL_MAX),
            safe(2.0)
        {
            if(eps <= 0)
            {
                Rf_error("eps must be positive in AccurateNumericDerivative");
            }
        }

        // for use as a function
        double
        operator()(double x) const;

    private:
        Fun& function;
        const double h;

        // sets maximum size of tableau
        const int ntab;

        // stepsize is decreased by factor con at each iteration
        const double con;
        const double con2;
        const double big;

        // return when error is safe worse than the best so far
        const double safe;
    };


// ***************************************************************************************************//

// get inverse second derivative of a function
template<class Fun>
    class NumericInvHessian
    {
    public:
        // save the function (object) reference and the epsilon parameter
        NumericInvHessian(Fun& function,
                          double eps = EPS) :
            function(function),
            eps(eps)
        {
        }

        // for use as a function
        double
        operator()(double x) const;

    private:
        Fun& function;
        const double eps;
    };

// ***************************************************************************************************//

// get more accurate inverse second derivative of a function
template<class Fun>
    class AccurateNumericInvHessian
    {
    public:
        // save the function (object) reference and the epsilon parameter
        AccurateNumericInvHessian(Fun& function) :
            function(function)
        {
        }

        // for use as a function
        double
        operator()(double x) const;

    private:
        Fun& function;
    };
// ***************************************************************************************************//

// for use as a function
template <class Fun>
double
CachedFunction<Fun>::operator()(double x)
{
    double ret = cache.getValue(x);

    // if not found in the old arguments
    if(R_IsNA(ret))
    {
        // we have to compute and save it.
        ret = function(x);
        cache.save(x, ret);
    }

    // finally return the value (either found or newly computed)
    return ret;
}

// ***************************************************************************************************//

// for use as a function
template <class Fun>
double
NumericDerivative<Fun>::operator()(double x) const
{
    // if abs(x) is (almost) zero, take delta = eps; otherwise delta = x*eps.
    // (important here: use of fabs instead of abs, which is an integer [!] version)
    double delta = (fabs(x) == 0.0) ? EPS : x * EPS;

    // return the estimated derivative at x:
    return (function(x + delta) - function(x)) / delta;
}

// ***************************************************************************************************//

// for use as a function
template <class Fun>
double
AccurateNumericDerivative<Fun>::operator()(double x) const
{
    // modified from
    // http://www.koders.com/cpp/fid192DE24D3250B7B3632C6F54490A67268C992E28.aspx?s=Chebyshev
    // (Numerical Recipes)

    double hh = h;
    double err = big;

    std::vector<MyDoubleVector> a (ntab, MyDoubleVector (ntab, 0));

    a[0][0] = (function(x + hh) - function(x - hh)) / (2.0 * hh);

    double answer = a[0][0];

    for (int i = 1; i < ntab; i++)
    {
        // Successive columns in the Neville tableau will go to smaller
        // stepsizes and higher orders of extrapolation
        hh /= con;

        // Try new smaller stepsize:
        a[0][i] = (function(x + hh) - function(x - hh)) / (2.0 * hh);

        double fac = con2;

        // Compute extrapolation of various orders, requiring no new function evaluations
        for (int j = 1; j <= i; j++)
        {
            a[j][i] = (a[j-1][i] * fac - a[j-1][i-1]) / (fac - 1.0);

            fac = con2 * fac;

            double errt = fmax(fabs(a[j][i] - a[j-1][i]),
                               fabs(a[j][i] - a[j-1][i-1]));

            // The error strategy is to compare each new extrapolation to one order lower,
            // both at the present stepsize and the previous one.

            // If error is decreased, save the improved answer.
            if (errt <= err)
            {
                err = errt;
                answer = a[j][i];
            }
        }

        // If higher order is worse by a significant factor safe, then quit early.
        if (fabs(a[i][i] - a[i-1][i-1]) >= safe * err)
            break;
    }

    return answer;
}

// ***************************************************************************************************//

// for use as a function
template <class Fun>
double
NumericInvHessian<Fun>::operator()(double x) const
{

    // if abs(x) is (almost) zero, take delta = eps; otherwise delta = x * eps.
    // (important here: use of fabs instead of abs, which is an integer [!] version)
    double delta = (fabs(x) == 0.0) ? eps : x * eps;
    delta = fabs(delta);

    // ensure that the difference between x and x + delta is exactly delta:
    double volatile temp = x + delta;
    delta = temp - x;

    // return the estimated inverse hessian at x:
    // (this is the analytical simplification of the centered differences formula, using the same
    // delta for the first and the second numerical derivative)
    // cf. Abramowitz/Stegun p. 884, formula 25.3.24.
    double forward = function(x + delta) - function(x);
    double backward = function(x) - function(x - delta);
    double ret = (delta * delta) / (forward - backward);

    return ret;
}

// ***************************************************************************************************//

// for use as a function
template <class Fun>
double
AccurateNumericInvHessian<Fun>::operator()(double x) const
{
    // get the accurate derivative of the function
    AccurateNumericDerivative<Fun> firstDerivative(function);

    // and the accurate derivative of the first derivative
    AccurateNumericDerivative<AccurateNumericDerivative<Fun> > secondDerivative(firstDerivative);

    // so the result is
    double ret = 1.0 / secondDerivative(x);

    // return that
    return ret;
}

// ***************************************************************************************************//

#endif /* FUNCTIONWRAPS_H_ */
