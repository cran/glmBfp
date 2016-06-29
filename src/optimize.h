/*
 * optimize.h
 *
 *  Created on: 08.12.2009
 *      Author: daniel
 */

#ifndef OPTIMIZE_H_
#define OPTIMIZE_H_

#include <types.h>

// modified from R's scalar function "optimize" routine "Brent_fmin"
// in /src/appl/fmin.c,
// in order to take a functor instead of an ordinary function, and does not
// take an info argument

// original documentation:

//    returns an approximation  x  to the point where  f  attains a minimum  on
//    the interval  (ax,bx)  is determined.
//
//    INPUT..
//
//    ax    left endpoint of initial interval
//    bx    right endpoint of initial interval
//    f     function which evaluates  f(x, info)  for any  x
//          in the interval  (ax,bx)
//    tol   desired length of the interval of uncertainty of the final
//          result ( >= 0.)
//
//    OUTPUT..
//
//    fmin  abcissa approximating the point where  f  attains a minimum
//
//        The method used is a combination of  golden  section  search  and
//    successive parabolic interpolation.  convergence is never much slower
//    than  that  for  a  Fibonacci search.  If  f  has a continuous second
//    derivative which is positive at the minimum (which is not  at  ax  or
//    bx),  then  convergence  is  superlinear, and usually of the order of
//    about  1.324....
//        The function  f  is never evaluated at two points closer together
//    than  eps*abs(fmin)+(tol/3), where eps is  approximately  the  square
//    root  of  the  relative  machine  precision.   if   f   is a unimodal
//    function and the computed values of   f   are  always  unimodal  when
//    separated  by  at least  eps*abs(x)+(tol/3), then  fmin  approximates
//    the abcissa of the global minimum of  f  on the interval  ax,bx  with
//    an error less than  3*eps*abs(fmin)+tol.  if   f   is  not  unimodal,
//    then fmin may approximate a local, but perhaps non-global, minimum to
//    the same accuracy.
//        This function subprogram is a slightly modified  version  of  the
//    Algol  60 procedure  localmin  given in Richard Brent, Algorithms for
//    Minimization without Derivatives, Prentice-Hall, Inc. (1973).


template <class Fun>
class Brent {
public:
    // setup the object
    // - function is the function to be maximized in the interval (lowerBound, upperBound)
    // - lowerBound: lower bound on the function argument
    // - upperBound: upper bound on the function argument
    // - precision: desired length of the interval of uncertainty of the final result ( > 0.0,
    //              default: fourth root of machine precision)
    Brent(Fun& function,
          double lowerBound,
          double upperBound,
          double precision=sqrt(EPS)) : // EPS is the square root of the machine precision
              function(function), lowerBound(lowerBound),
                    upperBound(upperBound), precision(precision)
        {
            // check arguments
            if (R_finite(lowerBound) == FALSE || R_finite(upperBound) == FALSE)
                Rf_error("Brent: bounds must be finite");

            if (lowerBound >= upperBound)
                Rf_error("Brent: lowerBound not smaller than upperBound");

            if (precision <= 0.0)
                Rf_error("Brent: precision not positive");
        }

    // minimize the function
    double
    minimize ();

private:
    // the function we want to maximize
    Fun& function;

    // constraints on the function argument
    const double lowerBound;
    const double upperBound;

    // settings
    const double precision;
};

// the actual minimization function
template<class Fun>
    double
    Brent<Fun>::minimize()
    {
        /*  c is the squared inverse of the golden ratio */
        const double c = (3.0 - sqrt(5.0)) * 0.5;

        /* Local variables */
        double a, b, d, e, p, q, r, u, v, w, x;
        double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

        /*  eps is approximately the square root of the relative machine precision. */
        eps = DOUBLE_EPS;
        tol1 = eps + 1.0; /* the smallest 1.000... > 1 */
        eps = sqrt(eps);

        a = lowerBound;
        b = upperBound;

        v = a + c * (b - a);
        w = v;
        x = v;

        d = 0.0; /* -Wall */
        e = 0.0;
        fx = function(x);
        fv = fx;
        fw = fx;
        tol3 = precision / 3.0;

        /*  main loop starts here ----------------------------------- */

        for (;;)
        {
            xm = (a + b) * 0.5;
            tol1 = eps * fabs(x) + tol3;
            t2 = tol1 * 2.0;

            /* check stopping criterion */

            if (fabs(x - xm) <= t2 - (b - a) * 0.5)
                break;
            p = 0.0;
            q = 0.0;
            r = 0.0;
            if (fabs(e) > tol1)
            { /* fit parabola */

                r = (x - w) * (fx - fv);
                q = (x - v) * (fx - fw);
                p = (x - v) * q - (x - w) * r;
                q = (q - r) * 2.0;
                if (q > 0.0)
                    p = - p;
                else
                    q = - q;
                r = e;
                e = d;
            }

            if (fabs(p) >= fabs(q * 0.5 * r) || p <= q * (a - x) || p >= q * (b - x))
            { /* a golden-section step */

                if (x < xm)
                    e = b - x;
                else
                    e = a - x;
                d = c * e;
            }
            else
            { /* a parabolic-interpolation step */

                d = p / q;
                u = x + d;

                /* f must not be evaluated too close to ax or bx */

                if (u - a < t2 || b - u < t2)
                {
                    d = tol1;
                    if (x >= xm)
                        d = - d;
                }
            }

            /* f must not be evaluated too close to x */

            if (fabs(d) >= tol1)
                u = x + d;
            else if (d > 0.)
                u = x + tol1;
            else
                u = x - tol1;

            fu = function(u);

            /*  update  a, b, v, w, and x */

            if (fu <= fx)
            {
                if (u < x)
                    b = x;
                else
                    a = x;
                v = w;
                w = x;
                x = u;
                fv = fw;
                fw = fx;
                fx = fu;
            }
            else
            {
                if (u < x)
                    a = u;
                else
                    b = u;
                if (fu <= fw || w == x)
                {
                    v = w;
                    fv = fw;
                    w = u;
                    fw = fu;
                }
                else if (fu <= fv || v == x || v == w)
                {
                    v = u;
                    fv = fu;
                }
            }
        }
        /* end of main loop */

        // return the found minimum
        return x;
    }



#endif /* OPTIMIZE_H_ */
