//
// Created by 马泽余 on 2019-11-01.
//

#include "powerSeries.h"


inline void clr(T *a, int l, int r)
{
    rep(i, l, r) a[i] = get(0, 0);
}

inline void cpy(T *a, T *b, int rr, int RR)
{
    rep(i, 0, rr) a[i] = b[i];
    rep(i, rr, RR) a[i] = get(0, 0);
}

inline void dft(T *a, int rr, bool inver = 0)
{
    int r = 0; while (1 << r < rr) r++;
    rep(i, 0, rr) if (i < rev[r][i]) swap(a[i], a[rev[r][i]]);
    for (int i = 0, ii = 1; ii < rr; i++, ii <<= 1)
    {
        T w = get(cos(2 * acos(-1) / (1 << (i + 1))), sin(2 * acos(-1) / (1 << (i + 1))));
        for (int j = 0; j < rr; j += (ii << 1))
        {
            T W = get(1, 0);
            for (int k = j; k < j + ii; k++)
            {
                T t = mul(a[k + ii], W);
                a[k + ii] = sub(a[k], t);
                a[k] = add(a[k], t);
                W = mul(W, w);
            }
        }
    }
    if (inver)
    {
        for (int i = 1; i < rr >> 1; i++) swap(a[i], a[rr - i]);
        rep(i, 0, rr) a[i] = div(a[i], rr);
    }
}


inline void mul(T *a, T *b, int rr, T *c)
{
    T *tmp1 = tmp[++cnt], *tmp2 = tmp[++cnt];
    cpy(tmp1, a, rr, rr << 1);
    dft(tmp1, rr << 1, 0);
    cpy(tmp2, b, rr, rr << 1);
    dft(tmp2, rr << 1, 0);
    rep(i, 0, rr << 1) c[i] = mul(tmp1[i], tmp2[i]);
    dft(c, rr << 1, 1);
    cnt -= 2;
}

T *multiply(T *a, T *b, int r) {
    int rr = 1 << r;
    rep (j, 0, r + 2) rep(i, 1, 1 << j) rev[j][i] = (rev[j][i >> 1] >> 1) | ((i & 1) << (j - 1));
    mul(a, b, rr, res);
    return res;
}


inline void mul2(T *a, T *b, int rr, T *c)
{
    T *tmp1 = tmp[++cnt];
    cpy(tmp1, a, rr, rr << 1);
    dft(tmp1, rr << 1, 0);
    rep(i, 0, rr << 1) c[i] = mul(tmp1[i], b[i]);
    dft(c, rr << 1, 1);
    cnt--;
}


inline void inverse_onestep(T *a, T *b, int rr)
{
    if (rr == 1)
    {
        a[0].r = 1 / b[0].r;
        a[0].i = 0;
        return;
    }
    T *tmp1 = tmp[++cnt];
    cpy(tmp1, b, rr, rr << 1);
    clr(a, rr >> 1, rr << 1);
    dft(a, rr << 1, 0);
    dft(tmp1, rr << 1, 0);
    rep(i, 0, rr << 1) a[i] = mul(a[i], sub(get(2, 0), mul(a[i], tmp1[i])));
    dft(a, rr << 1, 1);
    cnt--;
}

T *inverse(T *a, int r) {
    for (int i = 0; i <= r; i++) inverse_onestep(res, a, 1 << i);
    return res;
}

inline void derive(T *a, int rr)
{
    rep(i, 0, rr - 1) a[i] = mul(a[i + 1], get((i + 1), 0));
}

inline void integrate(T *a, int rr)
{
    for (int i = rr - 1; i >= 1; i--) a[i] = div(a[i - 1], i);
}

T *integ(T *a, int r) {
    cpy(res, a, 1 << r, 1 << r);
    integrate(res, 1 << r);
    return res;
}

T *derivative(T *a, int r) {
    cpy(res, a, 1 << r, 1 << r);
    derive(res, 1 << r);
    return res;
}


inline void ln(T *a, T *b, int rr, T *invb)
{
    T *tmp1 = tmp[++cnt];
    cpy(tmp1, b, rr, rr);
    derive(tmp1, rr);
    //inverse(a, b, rr);
    mul(invb, tmp1, rr, a);
    integrate(a, rr);
    a[0] = get(0, 0);
    cnt--;
}

T *ln(T *a, int r) {
    inverse(a, r);
    cpy(inv_xn, res, 1 << r, 1 << r);
    ln(res, a, 1 << r, inv_xn);
    return res;
}

inline void exp_onestep(T *a, T *b, int rr)
{
    if (rr == 1)
    {
        a[0] = get(1, 0);
        return;
    }
    T *tmp1 = tmp[++cnt], *tmp2 = tmp[++cnt];
    clr(a, rr >> 1, rr << 1);
    inverse_onestep(inv_xn_, a, rr >> 1);
    cpy(inv_xn, inv_xn_, rr >> 1, rr);
    inverse_onestep(inv_xn, a, rr);
    ln(tmp1, a, rr, inv_xn);
    clr(tmp1, rr, rr << 1);
    cpy(tmp2, b, rr, rr << 1);
    dft(a, rr << 1, 0);
    dft(tmp1, rr << 1, 0);
    dft(tmp2, rr << 1, 0);
    rep(i, 0, rr << 1) a[i] = mul(a[i], sub(add(get(1, 0), tmp2[i]), tmp1[i]));
    dft(a, rr << 1, 1);
    cnt -= 2;
}

T *exp(T *a, int r) {
    for (int i = 0; i <= r; i++) {
        exp_onestep(res, a, 1 << i);
    }
    return res;
}