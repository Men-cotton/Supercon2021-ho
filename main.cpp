#include <bits/stdc++.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
namespace sc21 {
#include "sc21.h"
}

// CとI_PROBを与えると誤差を返す関数。
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
            I[t + 1][i] = I[t][i] + sc21::BETA * S[t][i] * I[t][i] + sum -
                          sc21::GAMMA * I[t][i];
            R[t + 1][i] = R[t][i] + sc21::GAMMA * I[t][i];
        }
    }

    double loss = 0.0;
    for (int i = 0; i < N_GROUP; i++) {
        loss += (I[sc21::T][i] - I_PROB[i]) * (I[sc21::T][i] - I_PROB[i]);
    }

    return loss;
}

void backprop(int C[][N_GROUP], double I_PROB[], double C_back[][N_GROUP]) {
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
            I[t + 1][i] = I[t][i] + sc21::BETA * S[t][i] * I[t][i] + sum -
                          sc21::GAMMA * I[t][i];
            R[t + 1][i] = R[t][i] + sc21::GAMMA * I[t][i];
        }
    }

    double S_back[sc21::T + 1][N_GROUP] = {0},
                            I_back[sc21::T + 1][N_GROUP] = {0},
                            R_back[sc21::T + 1][N_GROUP] = {0};
    C_back[N_GROUP][N_GROUP] = {0};
    for (int i = 0; i < N_GROUP; i++) {
        I_back[sc21::T][i] = I[sc21::T][i] - I_PROB[i];
    }
    for (int t = sc21::T - 1; t >= 0; t--) {
        for (int i = 0; i < N_GROUP; i++) {
            R_back[t][i] = R_back[t + 1][i];
            I_back[t][i] = sc21::GAMMA * R_back[t + 1][i] + I_back[t + 1][i] +
                           sc21::BETA * I_back[t + 1][i] * S[t][i] -
                           sc21::GAMMA * I_back[t + 1][i] -
                           sc21::BETA * S[t][i] * S_back[t + 1][i];
            S_back[t][i] = sc21::BETA * I_back[t + 1][i] * I[t][i] +
                           S_back[t + 1][i] -
                           sc21::BETA * I[t][i] * S_back[t + 1][i];
            double sum = I_back[t + 1][i] - S_back[t + 1][i];
            S_back[t][i] += sc21::BETA2 * sum;
            sum *= sc21::BETA2 * S[t][i];

            for (int j = 0; j < N_GROUP; j++) {
                I_back[t][i] += sum * C[i][j];
                C_back[i][j] += sum * I[t][j];
            }
        }
    }
}
int main() {}