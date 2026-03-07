#include <Rcpp.h>
#include <vector>
#include <iostream>
using namespace Rcpp;
using namespace std;

NumericVector concatenate(const NumericVector &x, const NumericVector &y) {
  NumericVector result = x;
  int n = y.size();
  for (int i = 0; i < n; i++) {
    result.push_back(y[i]);
  }
  return result;
}

// [[Rcpp::export]]
List piecewise_quadratic(NumericMatrix coef, NumericVector b) {
    const NumericVector& A = coef.column(0);
    const NumericVector& B = coef.column(1);
    const NumericVector& C = coef.column(2);
    int m = b.length();
    double sum_A = 0, sum_B = 0, sum_C = 0;
    for (int i = 0; i < m; i++) {
        if (i == 0 || b[i-1] == R_PosInf) {
            sum_A += A[i];
            sum_B += B[i];
            sum_C += C[i];
        }
    }
    vector<pair<double, int> > end_pts;
    end_pts.reserve(m);
    for (int i = 0; i < m; i++) {
        end_pts.push_back(make_pair(b[i], i));
    }
    sort(end_pts.begin(), end_pts.end());
    double par = 0;
    double value = R_PosInf;
    for (int i = 0; i < m; i++) {
        double l = 0;
        if (i == 0) {
            l = R_NegInf;
        } else {
            l = end_pts[i-1].first;
        }
        double u = end_pts[i].first;
        if (l > R_NegInf) {
            double f = sum_A * l * l + sum_B * l + sum_C;
            if (f < value) {
                value = f;
                par = l;
            }
        }
        if (sum_A > 0 && -sum_B/2/sum_A > l && -sum_B/2/sum_A < u) {
            double f = sum_C - sum_B * sum_B / 4 / sum_A;
            if (f < value) {
                value = f;
                par = -sum_B / 2 / sum_A;
            }
        }
        if (u == R_PosInf) {
            break;   
        }
        int idx = end_pts[i].second;
        sum_A += A[idx+1] - A[idx];
        sum_B += B[idx+1] - B[idx];
        sum_C += C[idx+1] - C[idx];
    }
    return List::create(Named("par") = par, Named("value") = value);    
}

// [[Rcpp::export]]
List bracket_quadratic(NumericMatrix coef, NumericVector a = NumericVector::create(), NumericVector b = NumericVector::create()) {
  const NumericVector& A = coef.column(0);
  const NumericVector& B = coef.column(1);
  NumericVector C = coef.column(2);
  
  int m = A.length();

  if (a.length() == 0) {
    a.push_back(0);
  }
  
  if (a.length() == 1) {
    a = rep(a[0], m);
  }

  if (b.length() == 0) {
    b.push_back(R_NegInf);
    b.push_back(R_PosInf);
  }
  
  if (b.length() == 2) {
    b = rep_each(b, m);
  }
  
  double sum_a = 0;
  for (int i = 0; i < m; i++) {
    sum_a += a[i];
    C[i] -= a[i];
  }
  
  vector<pair<double, int> > end_pts;
  end_pts.reserve(m * 2);
  bool included[m];
  for (int i = 0; i < m; i++) {
    end_pts.push_back(make_pair(b[i], i));
    end_pts.push_back(make_pair(b[m+i], i));
    included[i] = false;
  }
  sort(end_pts.begin(), end_pts.end());
  double par = 0;
  double value = R_PosInf;
  double sum_A = 0, sum_B = 0, sum_C = 0;
  
  for (int i = 0; i < (int)end_pts.size() - 1; i++) {
    int p = end_pts[i].second;
    included[p] = !included[p];
    double sign = included[p] * 2 - 1;
    sum_A += sign * A[p];
    sum_B += sign * B[p];
    sum_C += sign * C[p];
    double l = end_pts[i].first;
    double u = end_pts[i+1].first;
    
    if (l < u) {
      if (l > R_NegInf) {
        double f = sum_A * l * l + sum_B * l + sum_C;
        if (f < value) {
          value = f;
          par = l;
        }
      }
      if (u < R_PosInf) {
        double f = sum_A * u * u + sum_B * u + sum_C;
        if (f < value) {
          value = f;
          par = u;
        }
      }
      if (sum_A > 0 && -sum_B/2/sum_A > l && -sum_B/2/sum_A < u) {
        double f = sum_C - sum_B * sum_B / 4 / sum_A;
        if (f < value) {
          value = f;
          par = -sum_B / 2 / sum_A;
        }
      }
    }
  }
  return List::create(Named("par") = par, Named("value") = value + sum_a);
}
