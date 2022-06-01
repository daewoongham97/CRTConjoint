#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List demean_center(NumericVector zz, int rows, int cols, NumericVector m) {
  NumericVector means(cols);
  for (int c = 0; c < cols; c++) {
    for (int r = 0; r < rows; r++) {
      int i = c* rows + r;
      zz[i] -= m[c];
    }
  }
  List out = List::create(_["vec"] = zz);
  return out;
}

// [[Rcpp::export]]
List fast_demean(NumericVector zz, int rows, int cols) {
  NumericVector means(cols);
  for (int c = 0; c < cols; c++) {
    double col_mean = 0.0;
    for (int r = 0; r < rows; r++) {
      int i = c* rows + r;
      col_mean += zz[i]/rows;
    }
    for (int r = 0; r < rows; r++) {
      int i = c* rows + r;
      zz[i] -= col_mean;
    }
    means[c] = col_mean;
  }
  List out = List::create(_["vec"] = zz, _["means"] = means);
  return out;
}




