/*
 * functionWraps.cpp
 *
 *  Created on: 06.11.2009
 *      Author: daniel
 */

#include <functionWraps.h>
#include <rcppExport.h>

using namespace Rcpp;

// ***************************************************************************************************//

// first R functions taking only a scalar

//// constructor
//RFunction::RFunction(SEXP R_fun)
//{
//    // allocate double scalar as function argument
//    Rf_protect(R_functionArgument = Rf_allocVector(REALSXP, 1));
//
//    // construct the function call
//    Rf_protect(R_functionCall = Rf_lcons(R_fun,
//                                  Rf_lcons(R_functionArgument,
//                                     R_NilValue)));
//
//    // get the function frame (is connected to the R function, so is already protected)
//    R_functionFrame = FRAME(R_fun);
//}

//// copy constructor
//RFunction::RFunction(const RFunction& old)
//{
//    // double scalar as function argument:
//    // take the same one as the old object.
//    Rf_protect(R_functionArgument = old.R_functionArgument);
//
//    // so also for the function call we must take the same one
//    Rf_protect(R_functionCall = old.R_functionCall);
//
//    // and then finally get the function frame
//    R_functionFrame = old.R_functionFrame;
//}


//// assignment operator
//RFunction& RFunction::operator=(const RFunction& rhs)
//{
//    if(&rhs != this)
//    {
//        // first unprotect the current members
//        Rf_unprotect(3);
//
//        // then assign the new members:
//
//        // double scalar as function argument:
//        // take the same one as the old object.
//        Rf_protect(R_functionArgument = rhs.R_functionArgument);
//
//        // so also for the function call we must take the same one
//        Rf_protect(R_functionCall = rhs.R_functionCall);
//
//        // and then finally get the function frame
//        R_functionFrame = rhs.R_functionFrame;
//    }
//
//    // return ourself
//    return *this;
//}
//

//// destructor
//RFunction::~RFunction()
//{
//    // unprotect the the function argument and the call
//    Rf_unprotect(2);
//
//    // everything else is automatically destroyed.
//}

// for use as a function
double
RFunction::operator()(double x) const
{
//    // input the argument
//    REAL(R_functionArgument)[0] = x;

//    // get the result
//    SEXP R_pars;
//    Rf_protect(R_pars = Rf_eval(R_functionCall,
//                                R_functionFrame));
//
//    // now convert to double vector
//    double ret = Rf_asReal(Rf_coerceVector(R_pars, REALSXP));
//
//    // unprotect local R result
//    Rf_unprotect(1);
//
//    // and return that
//    return ret;

    // shorter:
    // 8/11/2013: debug error
    // "'rho' must be an environment not pairlist: detected in C-level eval"
    // In gdb I find the following "classes" for the arguments:
    // (see Rinternals.h, line 98 ff.)

    double ret = as<double>(fun(x));
    return ret;

//    return Rf_asReal(Rf_coerceVector(Rf_eval(R_functionCall, // SEXP type 6: language constructs
//                                             R_functionFrame), // SEXP type 2: lists of dotted pairs
//                                     REALSXP));
}

// ***************************************************************************************************//


// now R functions taking a fixed length vector


//// constructor
//VectorRFunction::VectorRFunction(SEXP R_function, R_len_t argsize) :
//        argsize(argsize)
//{
//    // allocate vector as function argument
//    Rf_protect(R_functionArgument = Rf_allocVector(REALSXP, argsize));
//
//    // construct the function call
//    Rf_protect(R_functionCall = Rf_lcons(R_function,
//                                  Rf_lcons(R_functionArgument,
//                                     R_NilValue)));
//
//    // get the function frame (is connected to the R function, so is already protected)
//    R_functionFrame = FRAME(R_function);
//}
//
//
//// destructor
//VectorRFunction::~VectorRFunction()
//{
//    // unprotect the function argument and call
//    Rf_unprotect(2);
//
//    // everything else is automatically destroyed.
//}

// for use as a function
double
VectorRFunction::operator()(const double* vector) const
{
    // copy the vector
    NumericVector arg(argsize);
    std::copy(vector, vector + argsize, arg.begin());

    // and then evaluate the function call.
    double ret = as<double>(fun(arg));
    return ret;
}


// ***************************************************************************************************//

// save one pair of argument and function value
void
Cache::save(double arg, double val)
{
    args.push_back(arg);
    vals.push_back(val);
}

// clear the cache
void
Cache::clear()
{
    args.clear();
    vals.clear();
}

// query for a function value
double
Cache::getValue(double arg) const
{
    // search for the argument
    MyDoubleVector::const_iterator iterVals = vals.begin();
    for(MyDoubleVector::const_iterator
            iterArgs = args.begin();
            iterArgs != args.end();
            ++iterArgs, ++iterVals)
    {
        // match -> return corresponding value.
        if(*iterArgs == arg)
        {
            return *iterVals;
        }
    }

    // if we have not found the argument, return NA
    return R_NaReal;
}

// initialize from an R list
Cache::Cache(List& rcpp_list) :
        args(as<MyDoubleVector>(rcpp_list["args"])),
        vals(as<MyDoubleVector>(rcpp_list["vals"]))
{
    if(args.size() != vals.size())
    {
        Rf_error("Lengths of args and vals vectors in R list converted to Cache object not equal!");
    }
}

// convert to an Rcpp list
List
Cache::convert2list() const
{
    return List::create(_["args"] = getArgs(),
                        _["vals"] = getVals());
}


// ***************************************************************************************************//

// End of file.
