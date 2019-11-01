//
// Created by 马泽余 on 2019-11-01.
//

#ifndef POWERSERIESMANIPULATIONS_POWERSERIES_H
#define POWERSERIESMANIPULATIONS_POWERSERIES_H

#endif //POWERSERIESMANIPULATIONS_POWERSERIES_H

#include <stdio.h>
#include <math.h>
//#include <complex>

#define N 1 << 10
#define rep(i, l, r) for (int i = l; i < r; i++)
using namespace std;
typedef struct complex{
    double r, i;
} T;

T get(double a, double b) {
    T x;
    x.r = a;
    x.i = b;
    return x;
}

T mul(T a, T b) {
    T res;
    res.r = a.r * b.r - a.i * b.i;
    res.i = a.r * b.i + a.i * b.r;
    return res;
}

T div(T a, double k) {
    T res;
    res.r = a.r / k;
    res.i = a.i / k;
    return res;
}

T add(T a, T b) {
    T res;
    res.r = a.r + b.r;
    res.i = a.i + b.i;
    return res;
}

T sub(T a, T b) {
    T res;
    res.r = a.r - b.r;
    res.i = a.i - b.i;
    return res;
}

void swap(T &a, T&b) {
    double t;
    t = a.r;
    a.r = b.r;
    b.r = t;
    t = a.i;
    a.i = b.i;
    b.i = t;
}

int rev[21][N], cnt;
T a[N], b[N], exp_[N], inv_xn[N], inv_xn_[N], tmp[20][N];
T res[N];

T *multiply(T *a, T *b, int r);
T *inverse(T *a, int r);
T *integ(T *a, int r);
T *derivative(T *a, int r);
T *ln(T *a, int r);
T *exp(T *a, int r);