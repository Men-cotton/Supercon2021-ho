#include <bits/stdc++.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "sc21.h"

//CとI_PROBを与えると誤差を返す関数。
double simulator(int C_[][N_GROUP], double I_PROB_[]) {
    double S[T + 1][N_GROUP] = {0}, I[T + 1][N_GROUP] = {0},
                 R[T + 1][N_GROUP] = {0};
    for (int i = 0; i < N_GROUP; i++) {
        S[0][i] = N[i];
    }
    S[0][0] -= 1.0;
    I[0][0] = 1.0;

    for (int t = 0; t < T; t++) {
        for (int i = 0; i < N_GROUP; i++) {
            double sum = 0;
            for (int j = 0; j < N_GROUP; j++) {
                sum += C_[i][j] * I[t][j];
            }
            sum *= BETA2 * S[t][i];
            S[t + 1][i] = S[t][i] - BETA * S[t][i] * I[t][i] - sum;
            I[t + 1][i] =
                I[t][i] + BETA * S[t][i] * I[t][i] + sum - GAMMA * I[t][i];
            R[t + 1][i] = R[t][i] + GAMMA * I[t][i];
        }
    }

    double loss = 0.0;
    for (int i = 0; i < N_GROUP; i++) {
        loss += (I[T][i] - I_PROB_[i]) * (I[T][i] - I_PROB_[i]);
    }

    return loss;
}
int main() {
}