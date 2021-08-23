#include <bits/stdc++.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
namespace sc21 {
#include "sc21.h"
}

template <int L>
struct mat;

template <int L>
struct vec {
    double v[L];
    double operator[](const int i) const { return v[i]; }
    double& operator[](const int i) { return v[i]; }
    vec operator=(const vec& a) { memcpy(v, a.v, sizeof(a.v)); }
    vec operator=(double* a) {
        for (int i = 0; i < L; i++) {
            v[i] = a[i];
        }
    }
    vec& operator+=(const vec& a) {
        for (int i = 0; i < L; i++) {
            v[i] += a[i];
        }
    }
    vec& operator-=(const vec& a) {
        for (int i = 0; i < L; i++) {
            v[i] -= a[i];
        }
    }
    vec& operator*=(const vec& a) {
        for (int i = 0; i < L; i++) {
            v[i] *= a[i];
        }
    }
    vec& operator*=(const double a) {
        for (int i = 0; i < L; i++) {
            v[i] *= a;
        }
    }
    vec operator+(const vec& a) const {
        vec res = *this;
        res += a;
        return res;
    }
    vec operator-(const vec& a) const {
        vec res = *this;
        res -= a;
        return res;
    }
    vec operator*(const vec& a) const {
        vec res = *this;
        res *= a;
        return res;
    }
    vec operator*(const double a) const {
        vec res = *this;
        res *= a;
        return res;
    }
    vec mul(mat<L>& b);
    double dot(vec& b) {
        double res = 0.0;
        for (int i = 0; i < L; i++) {
            res += (*this)[i] * b[i];
        }
        return res;
    }
    void show() {
        for (int i = 0; i < L - 1; i++) {
            printf("%.5lf ", v[i]);
        }
        printf("%.5lf\n", v[L - 1]);
    }
};

template <int L>
struct mat {
    vec<L> m[L];
    vec<L> operator[](const int i) const { return m[i]; }
    vec<L>& operator[](const int i) { return m[i]; }
    mat operator=(const mat& a) { memcpy(m, a.m, sizeof(a.m)); }
    void show() {
        for (int i = 0; i < L; i++) {
            m[i].show();
        }
    }
};

//ベクトルと行列の積
template <int L>
vec<L> vec<L>::mul(mat<L>& b) {
    vec res = {};
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            res[i] += (*this)[j] * b[j][i];
        }
    }
    return res;
}

int main() {
    // sample 微妙
    double c[2] = {1.0, 2.0};
    vec<2> a, b = {1.0, 2.0};
    a = c;
    a *= 2.0;
    a.show();
    printf("%.5lf\n", a.dot(b));

    mat<2> m = {vec<2>{1.0, 2.0}, vec<2>{3.0, 4.0}};
    b = b.mul(m);
    b.show();
}