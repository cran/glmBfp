/*
 * coxfit.h
 *
 *  Created on: 03.12.2012
 *      Author: daniel
 */

#ifndef COXFIT_H_
#define COXFIT_H_

#include <types.h>

//  Dynamic_2d_array class by David Maisonave (609-345-1007) (www.axter.com)
//  Description:
//      The dynamic array class listed below is more efficient then other
//      similar classes that use temporary objects as return types, or use
//      an std::vector as a return type.
//
//      It's also more compatible with a C style 2D array, in that the
//      array is in one continuous memory block. This makes it possible
//      to pass this array object to a C Function that has a C-Style
//      2D array for a parameter.
//  Example usage:
/*

Dynamic_2d_array<int> MyIntArray(12, 34);
MyIntArray[0][1] = 123;
cout << MyIntArray[0][1] << endl;

*/

template < class T >
class Dynamic_2d_array
{
public:
  // constructor
  Dynamic_2d_array(int row, int col) :
    m_row(row),
    m_col(col),
    m_data((row!=0 && col!=0) ? new T[row*col] : NULL)
  {}

  // copy ctr
  Dynamic_2d_array(const Dynamic_2d_array& src) :
    m_row(src.m_row),
    m_col(src.m_col),
    m_data((src.m_row!=0 && src.m_col!=0) ? new T[src.m_row*src.m_col] : NULL)
  {
    for(int r=0; r<m_row; ++r)
      for(int c=0; c<m_col; ++c)
        (*this)[r][c] = src[r][c];
  }

  // destructor
  ~Dynamic_2d_array()
  {
    delete[] m_data;
  }

  // non-const access
  inline T* operator[](int i)
  {
    return (m_data + (m_col*i));
  }

  // const access
  inline T const*const operator[](int i) const
  {
    return (m_data + (m_col*i));
  }

private:
  const int m_row;
  const int m_col;
  T* m_data;
};

// Note that the class Dynamic_2d_array automatically allocates and
// deallocates the memory.
typedef Dynamic_2d_array<long> LongMatrix;
typedef Dynamic_2d_array<double> DoubleMatrix;
typedef Dynamic_2d_array<int> IntMatrix;

// **************************************************************************************


struct CoxfitResults
{
    // constructor:
    CoxfitResults(const int nCovs) :
        coefs(nCovs),
        imat(nCovs, nCovs),
        u(nCovs),
        loglik(2),
        flag(0)
    {
        // fill with initial values for the null model:
        coefs.fill(0.0);
    }


    // the coefficient vector
    AVector coefs;

    // the covariance matrix estimate
    AMatrix imat;

    // the score vector
    AVector u;

    // the log likelihoods of the initial coef vector and
    // the final one
    MyDoubleVector loglik;

    // convergence flag
    int flag;
};


// **************************************************************************************


// The class which does the Cox model fitting
class Coxfit {

public:

    // ctr
    Coxfit(const AVector& survTimes,
           const IntVector& censInd,
           const AMatrix& X,
           const AVector& weights,
           const AVector& offsets,
           const int method,
           double eps = 1e-09,
           double tolerChol = pow(DOUBLE_EPS, 0.75),
           int iterMax = 40,
           double tolerInf = 1e-05) :
               survTimes(survTimes),
               censInd(censInd),
               X(X),
               weights(weights),
               offsets(offsets),
               method(method),
               nObs(survTimes.size()),
               nCovs(X.n_cols),
               strata(nObs, 0),
               imat(nCovs, nCovs),
               results(nCovs),
               eps(eps),
               tolerChol(tolerChol),
               iterMax(iterMax),
               tolerInf(tolerInf)
               {}

    // getter for results
    CoxfitResults
    finalizeAndGetResults();

    // check results
    void
    checkResults() const;

    // the fit function, returns the number of
    // required iterations
    int
    fit();

    // compute the residual deviance of this model
    double
    computeResidualDeviance();

private:

    // inputs
    AVector survTimes;
    IntVector censInd;
    AMatrix X;
    AVector weights;
    AVector offsets;
    const int method;

    const int nObs;
    const int nCovs;

    IntVector strata;

    // temporary storage
    DoubleMatrix imat;

    // outputs
    CoxfitResults results;

    // numerical options
    const double eps;
    const double tolerChol;
    const int iterMax;
    const double tolerInf;
};

#endif /* COXFIT_H_ */
