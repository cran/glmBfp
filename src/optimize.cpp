/*
 * optimize.cpp
 *
 *  Created on: 08.12.2009
 *      Author: daniel
 */

#include <optimize.h>
#include <rcppExport.h>
#include <functionWraps.h>

using namespace Rcpp;

// ***************************************************************************************************//

// just an R interface to the brent optimization routine, for regression testing purposes.
// [[Rcpp::export]]
SEXP
cpp_optimize(SEXP R_function,  SEXP R_minx, SEXP R_maxx, SEXP R_precision)
{
    // extract function
    // R_interface = CDR(R_interface);
    // SEXP R_function = CAR(R_interface);
    // 
    // // constraints on the argument
    // R_interface = CDR(R_interface);
    // SEXP R_minx = CAR(R_interface);
    // 
    // R_interface = CDR(R_interface);
    // SEXP R_maxx= CAR(R_interface);
    // 
    // // and the settings
    // R_interface = CDR(R_interface);
    // SEXP R_precision = CAR(R_interface);

    // wrap the R function to a cached function
    RFunction fun(R_function);
    CachedFunction<RFunction> cachedFun(fun);

    // then minimize (with constraints)
    Brent<CachedFunction<RFunction> > brent(cachedFun,
                                            Rf_asReal(R_minx),
                                            Rf_asReal(R_maxx),
                                            Rf_asReal(R_precision));

    double xMin = brent.minimize();

    // now the inverse Hessian at the minimum
    AccurateNumericInvHessian<CachedFunction<RFunction> > funInvHess(cachedFun);
    double invHessMin = funInvHess(xMin);

    // pack results into R list
    return List::create(_["par"] = xMin,
                        _["inv.hessian"] = invHessMin,
                        _["evaluations"] = cachedFun.getCache().convert2list());
}

// ***************************************************************************************************//

// End of file.
