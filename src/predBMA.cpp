#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix predBMAcpp(NumericMatrix SurvMat, NumericMatrix LpMat, NumericVector WtVec) {
  
  int nModels = SurvMat.nrow();
  int nTimes = SurvMat.ncol();
  int nObs = LpMat.nrow();
  
  // arma::mat Pred(nTimes, nObs) //nTimes rows and nObs columns
  Rcpp::NumericMatrix Pred(nTimes, nObs); //nTimes rows and nObs columns
  
  for(int i=0; i < nModels; i++){
    for(int t=0; t < nTimes; t++){
      for(int Ob=0; Ob < nObs; Ob++){
        double S = SurvMat(i,t);
        double eLP = LpMat(Ob,i);
        Pred(t,Ob) = Pred(t,Ob) +  WtVec[i] * pow(S,eLP);
      }
    }
    //Rprintf("Model %d of %d",i,nModels);
  }
  
  return Pred;
}

