#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double constant_funC (NumericVector W, int K){
  double firstp=0;
  double secondp=0;
  double sumW=0;
  int n=W.size();

  for(int i=0; i<n; ++i){
    sumW +=W[i];
  }

  for(int i=1; i<=sumW; ++i){
    firstp +=log(i);
  }

  for(int k=1; k<=n; ++k){
    double total=0;
    if (W[k-1]==0){
      return total;
    }else{
      for(int i=1; i<=W[k-1]; ++i){
        total +=log(i);
      }
    }
    secondp +=total;
  }

  return firstp-secondp;
}
