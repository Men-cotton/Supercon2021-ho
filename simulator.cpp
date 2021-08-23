#include <bits/stdc++.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
namespace sc21{
    #include "sc21.h"
}

//CとI_PROBを与えると誤差を返す関数。
double simulator(int C[][N_GROUP], double I_PROB[]) {
    double S[sc21::T + 1][N_GROUP] = {0}, I[sc21::T + 1][N_GROUP] = {0},
                 R[sc21::T + 1][N_GROUP] = {0};
    for (int i = 0; i < N_GROUP; i++) {
        S[0][i] = sc21::N[i];
    }
    S[0][0] -= 1.0;
    I[0][0] = 1.0;

    for (int t = 0; t < sc21::T; t++) {
        for (int i = 0; i < N_GROUP; i++) {
            double sum = 0;
            for (int j = 0; j < N_GROUP; j++) {
                sum += C[i][j] * I[t][j];
            }
            sum *= sc21::BETA2 * S[t][i];
            S[t + 1][i] = S[t][i] - sc21::BETA * S[t][i] * I[t][i] - sum;
            I[t + 1][i] =
                I[t][i] + sc21::BETA * S[t][i] * I[t][i] + sum - sc21::GAMMA * I[t][i];
            R[t + 1][i] = R[t][i] + sc21::GAMMA * I[t][i];
        }
    }

    double loss = 0.0;
    for (int i = 0; i < N_GROUP; i++) {
        loss += (I[sc21::T][i] - I_PROB[i]) * (I[sc21::T][i] - I_PROB[i]);
    }

    return loss;
}
int main() {
}