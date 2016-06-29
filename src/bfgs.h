/*
 * bfgs.h
 *
 *  Created on: 13.11.2009
 *      Author: daniel
 */

#ifndef BFGS_H_
#define BFGS_H_

#include <functionWraps.h>

// ***************************************************************************************************//


// the main Bfgs class
template<class Fun>
class Bfgs {
public:
    // setup the object
    // - function is the function to be maximized (possibly constrained)
    // - verbose: print progress on R screen?
    // - lowerBound: lower bound on the function argument (default: -Inf)
    // - upperBound: upper bound on the function argument (default: +Inf)
    // - precision is precision parameter for the flatness at the maximum.
    // - ftol and gtol are convergence parameters
    // - stepsize is a parameter for the first ("crude") step of the algorithm
    Bfgs(Fun& function,
         bool verbose=false,
         double lowerBound=R_NegInf,
         double upperBound=R_PosInf,
         double precision=0.00001,
         double ftol=0.0001,
         double gtol=0.9,
         double stepsize=3.0) :
             function(function),
             functionDeriv(function),
             verbose(verbose),
             lowerBound(lowerBound),
             upperBound(upperBound),
             precision(precision),
             ftol(ftol),
             gtol(gtol),
             stepsize(stepsize)
             {
             }

    // minimize the function and compute the minimum and the inverted hessian at the minimum
    int
    minimize (double x0, double& xMin, double& invHessMin);

private:
    // the function we want to maximize
    Fun& function;

    // and its numerical derivative
    NumericDerivative<Fun> functionDeriv;

    // echo progress?
    const bool verbose;

    // constraints on the function argument
    const double lowerBound;
    const double upperBound;

    // settings
    const double precision;
    const double ftol;
    const double gtol;
    const double stepsize;


    // internal linesearch class:
    // (benefit is encapsulation and access to parent class members,
    // here function, functionDeriv, and settings)
    class Linesearch {
    public:
        Linesearch(const Bfgs& bfgs, // reference to instance of enclosing class must be passed!!
                   // there is no automatic connection of nested class instances to the enclosing class instance.
                   double x,
                   double direction) :
                       bfgs(bfgs),
                       x(x),
                       direction(direction),
                       maxAlpha((direction < 0) ?
                                ((bfgs.lowerBound - x) / direction) :
                                ((bfgs.upperBound - x) / direction))
                       {
                       }

        // do the linesearch for a given start point alpha1
        double
        operator()(double alpha1) const;

    private:

        // reference to the enclosing object
        const Bfgs& bfgs;

        // the point x and the direction of the search
        const double x;
        const double direction;

        // the maximum alpha to be chosen: set by the constructor.
        const double maxAlpha;

         // objective and derivative of objective function are member functions
        double
        phi(double alpha) const
        {
            return bfgs.function(x + alpha * direction);
        }

        double
        phi_(double alpha) const
        {
            return bfgs.functionDeriv(x + alpha * direction) * direction;
        }

        // other helper functions
        bool
        armijo(double alpha) const
        {
            return phi(alpha) <= (phi(0.0) + bfgs.ftol * alpha * phi_(0.0));
        }

        bool
        curvature(double alpha) const
        {
            return fabs(phi_(alpha)) <= - bfgs.gtol * phi_(0.0);
        }

        double
        interpolate(double x0,
                    double x1,
                    double f0,
                    double f1,
                    double g0) const;

        double
        zoom(double alpha_lo,
             double alpha_hi) const;
    };
};


// ***************************************************************************************************//

template <class Fun>
double
Bfgs<Fun>::Linesearch::zoom(double alpha_lo,
                       double alpha_hi) const
{
    for (int i = 0; i < 30; ++i)
    {
        double alpha = interpolate(alpha_lo, alpha_hi, phi(alpha_lo), phi(alpha_hi), phi_(alpha_lo));

        if ((! armijo(alpha)) || (phi(alpha) >= phi(alpha_lo)))
        {
            alpha_hi = alpha;
        }
        else
        {
            if (curvature(alpha))
            {
                return alpha;
            }

            // started going uphill?
            if (phi_(alpha) * (alpha_hi - alpha_lo) >= 0.0)
            {
                alpha_hi = alpha_lo;
            }
            alpha_lo = alpha;
        }
    }
    return 0.0; // not enough progress; give up.
}


// Implements the More-Thuente linesearch algorithm.  The notation is taken from
// Nocedal and Wright (2006).  This function returns an approximate solution to
//
//       argmin_{alpha \in [0, maxAlpha]} phi(alpha).
//
// To avoid getting stuck in local minima, the algorithm attempts to find a
// point that satisfies the Wolfe conditions:
//  - Armijo condition: the function has decreased enough.  This is controlled
// by the "ftol" parameter.
//  - curvature condition: the function has flattened out enough.  This is
// controlled by the "gtol" parameter.
// (It will succeed at this, unless machine precision gets in the way.)
//
// The basic idea of the algorithm is:
//  - do a crude search to find an interval that contains a local minimum
// that satisfies the Wolfe conditions.
//  - do a bisection-like search on that interval, until a point is found
// that satisfies the Wolfe conditions.
//
// In particular, the algorithm doesn't attempt to find a local minimizer.
// It just tries to find a rough guess, which is much faster.  The function
// that calls linesearch will do it many times.  Since the Wolfe conditions
// are relative to the starting point of "linesearch", each subsequent call
// will give tighter approximations.
//
// The algorithm will terminate if the objective eventually starts decreasing,
// or a finite constraint on alpha is specified.
//
// Parameters:
//  - phi : R -> R is the objective function to minimize.
//  - phi_ : R -> R is the derivative of phi.
//  - alpha1 is the starting point.
//  - maxAlpha is an upper-bound constraint on alpha (can be Inf).
template <class Fun>
double
Bfgs<Fun>::Linesearch::operator()(double alpha1) const
{
    // initial check
    if (! (phi_(0) < 0.0))
    {
        Rf_warning("\nBfgs: phi_(0) >= 0.0 in linesearch algorithm");
    }

    // initialize old and new alpha
    double alpha_ = 0.0;
    double alpha = (alpha1 >= maxAlpha) ? maxAlpha / 2.0 : alpha1;

    // iterate at most 100 times
    for (int i = 0; i < 100; ++i)
    {
        if ((i > 0 && (phi(alpha) >= phi(alpha_))) || (! armijo(alpha)))
        {
            return zoom(alpha_, alpha);
        }
        else if (curvature(alpha))
        {
            return alpha;
        }
        else if (phi_(alpha) >= 0.0)
        {
            return zoom(alpha, alpha_);
        }
        else
        {
            alpha_ = alpha;
            alpha = fmin((alpha + maxAlpha) / 2.0,
                         alpha * bfgs.stepsize);
        }
    }

    // return the last alpha if we did not break out of the loop beforehand.
    return alpha;
}

// ***************************************************************************************************//

// the actual minimization function:
// returns
// -1           lost precision
// 0            ok
// 1            change not large enough
template <class Fun>
int
Bfgs<Fun>::minimize(double x0, double& xMin, double& invHessMin)
{
    // first check that start value fulfills constraints.
    const bool insideBounds = (x0 >= lowerBound) && (x0 <= upperBound);
    if(! insideBounds)
    {
        Rf_error("Start value x0=%f for BFGS minimization not in admissible interval [%f, %f]", x0, lowerBound, upperBound);
    }

    // initialization
    xMin = x0;
    invHessMin = 1.0;

    int iter = 0;

    if(verbose)
    {
        Rprintf("\nBfgs: Starting BFGS minimization ...");
    }

    // as long as the gradient is far enough away from zero
    while (fabs(functionDeriv(xMin)) > precision)
    {
        ++iter;

        // possibly echo progress
        if (verbose)
        {
            Rprintf("\nBfgs: now at iteration %d", iter);
        }

        // minimize in the direction of p
        double p = - invHessMin * functionDeriv(xMin);

        // linesearch for factor alpha, starting from alpha = 1,
        // with point x into direction p.

        // first get maximum alpha


        double alpha = Linesearch(*this, xMin, p)(1.0);

        // check result of linesearch:
        // have we lost precision?
        if (alpha == 0.0)
        {
            // then warn
            if(verbose)
            {
                Rprintf("\nBfgs: Lost precision in linesearch of iteration %d", iter);
                Rprintf("\nBfgs: Finished minimization.");
            }

            // and give back the current results and the error code
            return(-1);
        }
        else if ((2.0 * fabs(alpha * p)) < precision)
        {
            // also break if the change is not large enough.
            if (verbose)
            {
                Rprintf("\nBfgs: Change not large enough in iteration %d", iter);
                Rprintf("\nBfgs: Finished minimization.");
            }

            return(1);
        }

        // compute the new x
        double x_ = xMin + alpha * p;

        // inverted Hessian is numeric deriv of g at the new x_:
        invHessMin = (x_ - xMin) / (functionDeriv(x_) - functionDeriv(xMin));

        // so now we are at x_
        xMin = x_;
    }

    if(verbose)
    {
        Rprintf("\nBfgs: Finished minimization.");
    }

    // return the convergence code
    return 0;
}

// ***************************************************************************************************//

// does quadratic interpolation to solve for the minimum on [x0, x1]
//       f(x) = ax^2 + bx + c
//
// (1) f(x0) = f0        => f0 = a x0^2 + b x0 + c
// (2) f(x1) = f1        => f1 = a x1^2 + b x1 + c
// (3) f'(x0) = g0       => g0 = 2a x0 + b
//
// (1)-(2)       f0 - f1 = a (x0^2 - x1^2) + b(x0 - x1)
// from (3):     b = g0 - 2a x0
// from (1)-(2): f0 - f1 = a (x0^2 - x1^2) + (g0 - 2a x0)(x0 - x1)
//               f0 - f1 - g0(x0 - x1) = a (x0^2 - x1^2) - 2a x0(x0 - x1)
//               f0 - f1 - g0(x0 - x1) = a [(x0^2 - x1^2) - 2 x0(x0 - x1)]
//               a = [f0 - f1 - g0(x0 - x1)] / [(x0^2 - x1^2) - 2 x0(x0 - x1)]
//                = - [f0 - f1 - g0(x0 - x1)] / [(x0 - x1)^2]
//
// Then, f'(x) = 2ax + b.  The "minimum" is at x = - b / (2a).
// (corrected here by DSB, it's not at x = -2a / b !!)
// of course, it is only the minimum if f''(x) = 2a > 0,
// which implies that f0 - f1 - g0(x0 - x1) < 0,
// or in other words g0 < (f0 - f1) / (x0 - x1) if x0 < x1.
template <class Fun>
double
Bfgs<Fun>::Linesearch::interpolate(double x0,
                                   double x1,
                                   double f0,
                                   double f1,
                                   double g0) const
{
    double dx = x0 - x1;
    double a = - (f0 - f1 - g0 * dx) / (dx * dx);
    double b = g0 - 2.0 * a * x0;

    double x = - b / (2.0 * a);

    // if it fails: return just the mid-point of [x0, x1]
    if ((R_finite(x) != 1) || (x < fmin(x0, x1)) || (x > fmax(x0, x1)))
    {
        return (x0 + x1) / 2.0;
    }
    else
    {
        return x;
    }
}

// ***************************************************************************************************//


#endif /* BFGS_H_ */
