/*
 * coxfit.cpp
 *
 *  Created on: 30.11.2012
 *      Author: daniel
 *
 *  Modified from R's internal cox fit function for Efron / Breslow ties method,
 *  R-2.15.2/src/library/Recommended/survival/src/coxfit6.c
 *  Three other functions related to Cholesky factorization and solving are
 *  also taken from the R sources.
 *
 */

#include <rcppExport.h>
#include <types.h>
#include <coxfit.h>

using namespace Rcpp;

/*
** A very few pathologic cases can cause the Newton Raphson iteration
**  path in coxph to generate a horrific argument to exp().  Since all these
**  calls to exp result in (essentially) relative risks we choose a
**  fixed value of LARGE on biological grounds: any number less than
**  1/(population of the earth) is essentially a zero, that is, an exponent
**  outside the range of +-23.
** A sensible numeric limit would be log(.Machine$double.xmax) which is
**  about 700, perhaps divided by 2 or log(n) to keep a few more bits.
**  However, passing this down the R calling chain to the c-routine is a lot
**  more hassle than I want to implement for this very rare case.
**
** Actually, the argument does not have to get large enough to have any
**  single exponential overflow.  In (start, stop] data we keep a running
**  sum of scores exp(x[i]*beta), which involves both adding subjects in and
**  subtracting them out.  An outlier x value that enters and then leaves can
**  erase all the digits of accuracy.  Most machines have about 16 digits of
**  accuracy and exp(21) uses up about 9 of them, leaving enough that the
**  routine doesn't fall on it's face.  (A user data set with outlier that
**  got exp(54) and a overlarge first beta on the first iteration led to this
**  paragraph.)  When beta-hat is infinite and x well behaved, the loglik
**  usually converges before xbeta gets to 15, so this protection should not
**  harm the iteration path of even edge cases; only fix those that truely
**  go astray.
**
** The truncation turns out not to be necessary for small values, since a risk
**  score of exp(-50) or exp(-1000) or 0 all give essentially the same effect.
** We only cut these off enough to avoid underflow.
*/

#define LARGE 22
#define SMALL -200

double coxsafe(double x) {
    if (x< SMALL) return(SMALL);
    if (x> LARGE) return(LARGE);
    return (x);
    }




/* $Id: cholesky2.c 11357 2009-09-04 15:22:46Z therneau $
**
** subroutine to do Cholesky decompostion on a matrix: C = FDF'
**   where F is lower triangular with 1's on the diagonal, and D is diagonal
**
** arguments are:
**     n         the size of the matrix to be factored
**     **matrix  a ragged array containing an n by n submatrix to be factored
**     toler     the threshold value for detecting "singularity"
**
**  The factorization is returned in the lower triangle, D occupies the
**    diagonal and the upper triangle is left undisturbed.
**    The lower triangle need not be filled in at the start.
**
**  Return value:  the rank of the matrix (non-negative definite), or -rank
**     it not SPD or NND
**
**  If a column is deemed to be redundant, then that diagonal is set to zero.
**
**   Terry Therneau
*/

int cholesky2(DoubleMatrix& matrix, int n, double toler)
    {
    double temp;
    int  i,j,k;
    double eps, pivot;
    int rank;
    int nonneg;

    nonneg=1;
    eps =0;
    for (i=0; i<n; i++) {
        if (matrix[i][i] > eps)  eps = matrix[i][i];
        for (j=(i+1); j<n; j++)  matrix[j][i] = matrix[i][j];
        }
    eps *= toler;

    rank =0;
    for (i=0; i<n; i++) {
        pivot = matrix[i][i];
        if (pivot < eps) {
            matrix[i][i] =0;
            if (pivot < -8*eps) nonneg= -1;
            }
        else  {
            rank++;
            for (j=(i+1); j<n; j++) {
                temp = matrix[j][i]/pivot;
                matrix[j][i] = temp;
                matrix[j][j] -= temp*temp*pivot;
                for (k=(j+1); k<n; k++) matrix[k][j] -= temp*matrix[k][i];
                }
            }
        }
    return(rank * nonneg);
    }

/* $Id: chinv2.c 11357 2009-09-04 15:22:46Z therneau $
**
** matrix inversion, given the FDF' cholesky decomposition
**
** input  **matrix, which contains the chol decomp of an n by n
**   matrix in its lower triangle.
**
** returned: the upper triangle + diagonal contain (FDF')^{-1}
**            below the diagonal will be F inverse
**
**  Terry Therneau
*/

void chinv2(DoubleMatrix& matrix , int n)
     {
     double temp;
     int i,j,k;

     /*
     ** invert the cholesky in the lower triangle
     **   take full advantage of the cholesky's diagonal of 1's
     */
     for (i=0; i<n; i++){
          if (matrix[i][i] >0) {
              matrix[i][i] = 1/matrix[i][i];   /*this line inverts D */
              for (j= (i+1); j<n; j++) {
                   matrix[j][i] = -matrix[j][i];
                   for (k=0; k<i; k++)     /*sweep operator */
                        matrix[j][k] += matrix[j][i]*matrix[i][k];
                   }
              }
          }

     /*
     ** lower triangle now contains inverse of cholesky
     ** calculate F'DF (inverse of cholesky decomp process) to get inverse
     **   of original matrix
     */
     for (i=0; i<n; i++) {
          if (matrix[i][i]==0) {  /* singular row */
                for (j=0; j<i; j++) matrix[j][i]=0;
                for (j=i; j<n; j++) matrix[i][j]=0;
                }
          else {
              for (j=(i+1); j<n; j++) {
                   temp = matrix[j][i]*matrix[j][j];
                   if (j!=i) matrix[i][j] = temp;
                   for (k=i; k<j; k++)
                        matrix[i][k] += temp*matrix[j][k];
                   }
              }
          }
     }

/*  $Id: chsolve2.c 11376 2009-12-14 22:53:57Z therneau $
**
** Solve the equation Ab = y, where the cholesky decomposition of A and y
**   are the inputs.
**
** Input  **matrix, which contains the chol decomp of an n by n
**   matrix in its lower triangle.
**        y[n] contains the right hand side
**
**  y is overwriten with b
**
**  Terry Therneau
*/

void chsolve2(const DoubleMatrix& matrix, int n, AVector& y)
     {
     int i,j;
     double temp;

     /*
     ** solve Fb =y
     */
     for (i=0; i<n; i++) {
          temp = y[i] ;
          for (j=0; j<i; j++)
               temp -= y[j] * matrix[i][j] ;
          y[i] = temp ;
          }
     /*
     ** solve DF'z =b
     */
     for (i=(n-1); i>=0; i--) {
          if (matrix[i][i]==0)  y[i] =0;
          else {
              temp = y[i]/matrix[i][i];
              for (j= i+1; j<n; j++)
                   temp -= y[j]*matrix[j][i];
              y[i] = temp;
              }
          }
     }




// getter for results
CoxfitResults
Coxfit::finalizeAndGetResults()
{
    // copy the contents of the DoubleMatrix imat into
    // the Armadillo matrix results.imat
    for(int i = 0; i < nCovs; ++i)
    {
        for(int j = 0; j < nCovs; ++j)
        {
            results.imat.at(i, j) = imat[i][j];
        }
    }
    // the other results are already up to date, because
    // their storage is directly used in the fit function.

    // then return the results
    return results;
}


// check results
void
Coxfit::checkResults() const
{
    // todo: replace errors by exceptions later on

    if(results.flag < nCovs)
    {
        Rf_error("Singular model!");
    }

    AVector infs = results.imat * results.u;

    if(iterMax > 1)
    {
        if (results.flag == 1000)
        {
            Rf_error("Ran out of iterations and did not converge");
        }
        else
        {
            IntVector notConverged;

            for(int i=0; i < nCovs; ++i)
            {
                if((infs[i] > eps) & (infs[i] > tolerInf * fabs(results.coefs[i])))
                {
                    notConverged.push_back(i + 1);
                }
            }

            if(notConverged.size() > 0)
            {
                Rf_warning("Loglik converged before some variables; beta may be infinite. ");
            }
        }
    }
}


// the fit function, returns the number of
// required iterations
int
Coxfit::fit()
{
    int i,j,k, person;

    double  wtave;
    double  denom=0, zbeta, risk;
    double  temp, temp2;
    int     ndead;  /* actually, the sum of their weights */
    double  newlk = 0.0;
    double  dtime, d2;
    double  deadwt;  /*sum of case weights for the deaths*/
    double  efronwt; /* sum of weighted risk scores for the deaths*/
    int     halving;    /*are we doing step halving at the moment? */
    int     nrisk;   /* number of subjects in the current risk set */

    /* copies of scalar input arguments */
    int     nused, nvar, maxiter;
    double  toler;
    int doscale;

    /* returned objects */
    int iter;

    /* get local copies of some input args */
    nused = nObs;
    nvar  = nCovs;
    maxiter = iterMax;
    toler = tolerChol;  /* tolerance for cholesky */
    doscale = 1;

    /*
     **  Set up the ragged arrays and scratch space
     */

    AVector a(nvar);
    AVector newbeta(nvar);
    AVector a2(nvar);
    AVector scale(nvar);

    DoubleMatrix cmat(nvar, nvar);
    DoubleMatrix cmat2(nvar, nvar);

    /*
     ** Subtract the mean from each covar, as this makes the regression
     **  much more stable.
     */
    for (i=0; i<nvar; i++) {
        temp=0;
        for (person=0; person<nused; person++) temp += X(person, i);
        temp /= nused;
        for (person=0; person<nused; person++)
        {
            X(person, i) -=temp;
        }
        if (doscale==1) {  /* and also scale it */
            temp =0;
            for (person=0; person<nused; person++) {
                temp += fabs(X(person, i));
            }
            if (temp > 0) temp = nused/temp;   /* scaling */
            else temp=1.0; /* rare case of a constant covariate */
            scale[i] = temp;
            for (person=0; person<nused; person++)
            {
                X(person, i) *= temp;
            }
        }
    }
    if (doscale==1) {
        for (i=0; i<nvar; i++) results.coefs[i] /= scale[i]; /*rescale initial betas */
    }
    else {
        for (i=0; i<nvar; i++) scale[i] = 1.0;
    }

    /*
     ** do the initial iteration step
     */
    strata[nused-1] =1;
    results.loglik[1] =0;
    for (i=0; i<nvar; i++) {
        results.u[i] =0;
        a2[i] =0;
        for (j=0; j<nvar; j++) {
            imat[i][j] =0 ;
            cmat2[i][j] =0;
        }
    }

    for (person=nused-1; person>=0; ) {
        if (strata[person] == 1) {
            nrisk =0 ;
            denom = 0;
            for (i=0; i<nvar; i++) {
                a[i] = 0;
                for (j=0; j<nvar; j++) cmat[i][j] = 0;
            }
        }

        dtime = survTimes[person];
        ndead =0; /*number of deaths at this time point */
        deadwt =0;  /* sum of weights for the deaths */
        efronwt=0;  /* sum of weighted risks for the deaths */
        while(person >=0 &&survTimes[person]==dtime) {
            /* walk through the this set of tied times */
            nrisk++;
            zbeta = offsets[person];    /* form the term beta*z (vector mult) */
            for (i=0; i<nvar; i++)
                zbeta += results.coefs[i]*X(person, i);
            zbeta = coxsafe(zbeta);
            risk = exp(zbeta) * weights[person];
            denom += risk;

            /* a is the vector of weighted sums of x, cmat sums of squares */
            for (i=0; i<nvar; i++) {
                a[i] += risk*X(person, i);
                for (j=0; j<=i; j++)
                    cmat[i][j] += risk * X(person, i) * X(person, j);
            }

            if (censInd[person]==1) {
                ndead++;
                deadwt += weights[person];
                efronwt += risk;
                results.loglik[1] += weights[person]*zbeta;

                for (i=0; i<nvar; i++)
                    results.u[i] += weights[person]*X(person, i);
                if (method==1) { /* Efron */
                    for (i=0; i<nvar; i++) {
                        a2[i] +=  risk*X(person, i);
                        for (j=0; j<=i; j++)
                            cmat2[i][j] += risk* X(person, i) * X(person, j);
                    }
                }
            }

            person--;
            if (strata[person]==1) break;  /*ties don't cross strata */
        }


        if (ndead >0) {  /* we need to add to the main terms */
            if (method==0) { /* Breslow */
                results.loglik[1] -= deadwt* log(denom);

                for (i=0; i<nvar; i++) {
                    temp2= a[i]/ denom;  /* mean */
                    results.u[i] -=  deadwt* temp2;
                    for (j=0; j<=i; j++)
                        imat[j][i] += deadwt*(cmat[i][j] - temp2*a[j])/denom;
                }
            }
            else { /* Efron */
                /*
                 ** If there are 3 deaths we have 3 terms: in the first the
                 **  three deaths are all in, in the second they are 2/3
                 **  in the sums, and in the last 1/3 in the sum.  Let k go
                 **  from 0 to (ndead -1), then we will sequentially use
                 **     denom - (k/ndead)*efronwt as the denominator
                 **     a - (k/ndead)*a2 as the "a" term
                 **     cmat - (k/ndead)*cmat2 as the "cmat" term
                 **  and reprise the equations just above.
                 */
                for (k=0; k<ndead; k++) {
                    temp = (double)k/ ndead;
                    wtave = deadwt/ndead;
                    d2 = denom - temp*efronwt;
                    results.loglik[1] -= wtave* log(d2);
                    for (i=0; i<nvar; i++) {
                        temp2 = (a[i] - temp*a2[i])/ d2;
                        results.u[i] -= wtave *temp2;
                        for (j=0; j<=i; j++)
                            imat[j][i] +=  (wtave/d2) *
                            ((cmat[i][j] - temp*cmat2[i][j]) -
                                    temp2*(a[j]-temp*a2[j]));
                    }
                }

                for (i=0; i<nvar; i++) {
                    a2[i]=0;
                    for (j=0; j<nvar; j++) cmat2[i][j]=0;
                }
            }
        }
    }   /* end  of accumulation loop */
    results.loglik[0] = results.loglik[1]; /* save the loglik for iter 0 */

    /* am I done?
     **   update the betas and test for convergence
     */
    for (i=0; i<nvar; i++) /*use 'a' as a temp to save u0, for the score test*/
        a[i] = results.u[i];

    results.flag= cholesky2(imat, nvar, toler);
    chsolve2(imat,nvar,a);        /* a replaced by  a *inverse(i) */


    /*
     **  Never, never complain about convergence on the first step.  That way,
     **  if someone HAS to they can force one iter at a time.
     */
    for (i=0; i<nvar; i++) {
        newbeta[i] = results.coefs[i] + a[i];
    }
    if (maxiter==0) {
        chinv2(imat,nvar);
        for (i=0; i<nvar; i++) {
            results.coefs[i] *= scale[i];  /*return to original scale */
            results.u[i] /= scale[i];
            imat[i][i] *= scale[i]*scale[i];
            for (j=0; j<i; j++) {
                imat[j][i] *= scale[i]*scale[j];
                imat[i][j] = imat[j][i];
            }
        }
        return 0;
    }

    /*
     ** here is the main loop
     */
    halving =0 ;             /* =1 when in the midst of "step halving" */
    for (iter=1; iter<= maxiter; iter++) {
        newlk =0;
        for (i=0; i<nvar; i++) {
            results.u[i] =0;
            for (j=0; j<nvar; j++)
                imat[i][j] =0;
        }

        /*
         ** The data is sorted from smallest time to largest
         ** Start at the largest time, accumulating the risk set 1 by 1
         */
        for (person=nused-1; person>=0; ) {
            if (strata[person] == 1) { /* rezero temps for each strata */
                denom = 0;
                nrisk =0;
                for (i=0; i<nvar; i++) {
                    a[i] = 0;
                    for (j=0; j<nvar; j++) cmat[i][j] = 0;
                }
            }

            dtime = survTimes[person];
            deadwt =0;
            ndead =0;
            efronwt =0;
            while(person>=0 && survTimes[person]==dtime) {
                nrisk++;
                zbeta = offsets[person];
                for (i=0; i<nvar; i++)
                    zbeta += newbeta[i]*X(person, i);
                zbeta = coxsafe(zbeta);
                risk = exp(zbeta) * weights[person];
                denom += risk;

                for (i=0; i<nvar; i++) {
                    a[i] += risk*X(person, i);
                    for (j=0; j<=i; j++)
                        cmat[i][j] += risk*X(person, i)*X(person, j);
                }

                if (censInd[person]==1) {
                    ndead++;
                    deadwt += weights[person];
                    newlk += weights[person] *zbeta;
                    for (i=0; i<nvar; i++)
                        results.u[i] += weights[person] *X(person, i);
                    if (method==1) { /* Efron */
                        efronwt += risk;
                    for (i=0; i<nvar; i++) {
                        a2[i] +=  risk*X(person, i);
                        for (j=0; j<=i; j++)
                            cmat2[i][j] += risk*X(person, i)*X(person, j);
                    }
                    }
                }

                person--;
                if (strata[person]==1) break; /*tied times don't cross strata*/
            }

            if (ndead >0) {  /* add up terms*/
                if (method==0) { /* Breslow */
                    newlk -= deadwt* log(denom);
                    for (i=0; i<nvar; i++) {
                        temp2= a[i]/ denom;  /* mean */
                        results.u[i] -= deadwt* temp2;
                        for (j=0; j<=i; j++)
                            imat[j][i] +=  (deadwt/denom)*
                            (cmat[i][j] - temp2*a[j]);
                    }
                }
                else  { /* Efron */
                    for (k=0; k<ndead; k++) {
                        temp = (double)k / ndead;
                        wtave= deadwt/ ndead;
                        d2= denom - temp* efronwt;
                        newlk -= wtave* log(d2);
                        for (i=0; i<nvar; i++) {
                            temp2 = (a[i] - temp*a2[i])/ d2;
                            results.u[i] -= wtave*temp2;
                            for (j=0; j<=i; j++)
                                imat[j][i] +=  (wtave/d2)*
                                ((cmat[i][j] - temp*cmat2[i][j]) -
                                        temp2*(a[j]-temp*a2[j]));
                        }
                    }

                    for (i=0; i<nvar; i++) { /*in anticipation */
                        a2[i] =0;
                        for (j=0; j<nvar; j++) cmat2[i][j] =0;
                    }
                }
            }
        }   /* end  of accumulation loop  */

        /* am I done?
         **   update the betas and test for convergence
         */
        results.flag = cholesky2(imat, nvar, toler);

        if (fabs(1-(results.loglik[1]/newlk))<= eps && halving==0) { /* all done */
            results.loglik[1] = newlk;
            chinv2(imat, nvar);     /* invert the information matrix */
            for (i=0; i<nvar; i++) {
                results.coefs[i] = newbeta[i]*scale[i];
                results.u[i] /= scale[i];
                imat[i][i] *= scale[i]*scale[i];
                for (j=0; j<i; j++) {
                    imat[j][i] *= scale[i]*scale[j];
                    imat[i][j] = imat[j][i];
                }
            }
            return iter;
        }

        if (iter== maxiter) break;  /*skip the step halving calc*/

        if (newlk < results.loglik[1])   {    /*it is not converging ! */
            halving =1;
            for (i=0; i<nvar; i++)
                newbeta[i] = (newbeta[i] + results.coefs[i]) /2; /*half of old increment */
        }
        else {
            halving=0;
            results.loglik[1] = newlk;
            chsolve2(imat,nvar,results.u);
            j=0;
            for (i=0; i<nvar; i++) {
                results.coefs[i] = newbeta[i];
                newbeta[i] = newbeta[i] + results.u[i];
            }
        }
    }   /* return for another iteration */

    /*
     ** We end up here only if we ran out of iterations
     */
    results.loglik[1] = newlk;
    chinv2(imat, nvar);
    for (i=0; i<nvar; i++) {
        results.coefs[i] = newbeta[i]*scale[i];
        results.u[i] /= scale[i];
        imat[i][i] *= scale[i]*scale[i];
        for (j=0; j<i; j++) {
            imat[j][i] *= scale[i]*scale[j];
            imat[i][j] = imat[j][i];
        }
    }
    results.flag = 1000;

    return iter;
}


// ***************************************************************************************************//

// compute the residual deviance of this model
double
Coxfit::computeResidualDeviance()
{
    // do the fitting
    fit();

    // get results: need to do that first!
    CoxfitResults fit = finalizeAndGetResults();

    // check results
    checkResults();

    // return residual deviance
    double ret = - 2.0 * (fit.loglik[0] - fit.loglik[1]);
    return ret;
}

// ***************************************************************************************************//

// 03/07/2013: add offsets

// just an R interface to the coxfit routine, for regression testing purposes.
// [[Rcpp::export]]
SEXP
cpp_coxfit(SEXP R_survTimes, SEXP R_censInd, SEXP R_offsets, SEXP R_X, SEXP R_method)
{
    // ---------------
    // get R objects:

    // // extract survival times
    // R_interface = CDR(R_interface);
    // SEXP R_survTimes = CAR(R_interface);
    // 
    // // censoring status
    // R_interface = CDR(R_interface);
    // SEXP R_censInd = CAR(R_interface);
    // 
    // // offsets
    // R_interface = CDR(R_interface);
    // SEXP R_offsets = CAR(R_interface);
    // 
    // // design matrix
    // R_interface = CDR(R_interface);
    // SEXP R_X= CAR(R_interface);
    // 
    // // and the tie method
    // R_interface = CDR(R_interface);
    // SEXP R_method = CAR(R_interface);

    // ---------------
    // unpack R objects:

    // survival times
    const NumericVector n_survTimes = R_survTimes;
    //const AVector survTimes(n_survTimes.begin(), n_survTimes.size(),
    //                      false);
    
    // errors with const and stuff... what if we copy into new memory?
    const AVector survTimes(n_survTimes.begin(), n_survTimes.size());

    // censoring status
    const IntVector censInd = as<IntVector>(R_censInd);

    // offsets
    const AVector offsets = as<NumericVector>(R_offsets);

    // design matrix
    const NumericMatrix n_X = R_X;
    //const AMatrix X(n_X.begin(), n_X.nrow(),
    //                n_X.ncol(), false);
    
    //Same issue as above L717:721
    const AMatrix X(n_X.begin(), n_X.nrow(),
                    n_X.ncol());
    

    // tie method
    const int method = as<int>(R_method);

    // ---------------
    // assign remaining arguments for Coxfit

    const int nObs = survTimes.size();

    const AVector weights = arma::ones<AVector>(nObs);

    // ---------------

    // get new Coxfit object
    Coxfit cox(survTimes,
               censInd,
               X,
               weights,
               offsets,
               method);

    // do the fitting
    const int nIter = cox.fit();

    // get results: need to do that first!
    CoxfitResults fit = cox.finalizeAndGetResults();

    // check results
    cox.checkResults();

    // pack results into R list and return that
    return List::create(_["coef"] = fit.coefs,
                        _["imat"] = fit.imat,
                        _["loglik"] = fit.loglik,
                        _["nIter"] = nIter);

}

