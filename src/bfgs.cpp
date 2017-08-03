/*
 * bfgs.cpp
 *
 *  Created on: 11.11.2009
 *      Author: daniel
 */

#include <bfgs.h>
#include <rcppExport.h>
#include <functionWraps.h>

using namespace Rcpp;

// ***************************************************************************************************//



// just an R interface to bfgs, for regression testing purposes.
// [[Rcpp::export]]
SEXP
cpp_bfgs(SEXP r_startval, SEXP r_function, SEXP r_minx, SEXP r_maxx, SEXP r_precision, SEXP r_verbose  )
{
    // extract function and start value for minimization
    // r_interface = CDR(r_interface);
    // SEXP r_startval = CAR(r_interface);
    // 
    // r_interface = CDR(r_interface);
    // SEXP r_function = CAR(r_interface);
    // 
    // // constraints on the argument
    // r_interface = CDR(r_interface);
    // SEXP r_minx = CAR(r_interface);
    // 
    // r_interface = CDR(r_interface);
    // SEXP r_maxx= CAR(r_interface);
    // 
    // // and the settings
    // r_interface = CDR(r_interface);
    // SEXP r_precision = CAR(r_interface);
    // 
    // r_interface = CDR(r_interface);
    // SEXP r_verbose = CAR(r_interface);

    // wrap the R function to a cached function
    RFunction fun(r_function);
    CachedFunction<RFunction> cachedFun(fun);

    // then minimize (with constraints)
    Bfgs<CachedFunction<RFunction> > bfgs(cachedFun,
                                          Rf_asLogical(r_verbose) == 1,
                                          Rf_asReal(r_minx),
                                          Rf_asReal(r_maxx),
                                          Rf_asReal(r_precision));

    double xMin = 0.0;
    double invHessMin = 0.0;

    int code = bfgs.minimize(Rf_asReal(r_startval),
                             xMin,
                             invHessMin);

    // pack results into R list
    return List::create(_["par"] = xMin,
                        _["inv.hessian"] = invHessMin,
                        _["evaluations"] = cachedFun.getCache().convert2list(),
                        _["code"] = code);
}

// ***************************************************************************************************//

// End of file.

