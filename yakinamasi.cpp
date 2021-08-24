#include <bits/stdc++.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

namespace sc21 {
#include "sc21.h"
}

unsigned int xor128() {
    static unsigned int x = 123456789, y = 362436069, z = 521288629,
                        w = 88675123;
    unsigned int t;
    t = (x ^ (x << 11));
    x = y;
    y = z;
    z = w;
    return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
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

void modify(int C[][N_GROUP],int change){
    if(change<=0)return;
    for(int i=0;i<change;i++){
        int a=xor128()%N_GROUP,b=xor128()%N_GROUP,c=xor128()%N_GROUP,d=xor128()%N_GROUP;
        if(a<b&&c<d&&C[a][b]==1&&C[c][d]==0){
            C[a][b]=0;
            C[b][a]=0;
            C[c][d]=1;
            C[d][c]=1;
            break;
        }
    }
}

const unsigned long long int cycle_per_sec = 2800000000;
 
unsigned long long int getCycle() {
    unsigned int low, high;
    __asm__ volatile ("rdtsc" : "=a" (low), "=d" (high));
    return ((unsigned long long int)low) | ((unsigned long long int)high << 32);
}
 
double getTime (unsigned long long int begin_cycle) {
    return (double)(getCycle() - begin_cycle) / cycle_per_sec;
}

// 焼きなまし法
void sa(int C[][N_GROUP]){
    auto start_cycle = getCycle();

    //初期化
    for(int i=0;i<N_GROUP;i++){
        for(int j=0;j<N_GROUP;j++){
            C[i][j]=0;
        }
    }
    for(int i=0;i<N_GROUP;i++){
        C[i][(i+1)%N_GROUP]=1;
        C[(i+1)%N_GROUP][i]=1;
    }
    for(int i=0;i<sc21::N_LINK-N_GROUP;i++){
        while(true){
            int a=xor128()%N_GROUP,b=xor128()%N_GROUP;
            if(a<b&&C[a][b]==0){
                C[a][b]=1;
                C[b][a]=1;
                break;
            }
        }
    }

    double start_temp=100000,end_temp=0;
    double TIME_LIMIT=9.50;

    while(true){
        double now_t = getTime(start_cycle);
        if(now_t>TIME_LIMIT)break;

        int new_state[N_GROUP][N_GROUP];

        for(int i=0;i<N_GROUP;i++){
            for(int j=0;j<N_GROUP;j++){
                new_state[i][j]=C[i][j];
            }
        }

        //だんだん変更量を減少させるように
        int change=1+10-10*(now_t)/TIME_LIMIT;
        
        modify(new_state,change);

        double new_score=simulator(new_state,sc21::I_PROB);
        double pre_score=simulator(C,sc21::I_PROB);

        double temp=start_temp+(end_temp-start_temp)*(sc21::TIME0-start_cycle)/TIME_LIMIT;
        double prob=exp((new_score-pre_score)/temp);

        if(prob>(xor128()%1280000)/1280000){
            for(int i=0;i<N_GROUP;i++){
                for(int j=0;j<N_GROUP;j++){
                    C[i][j]=new_state[i][j];
                }
            }
        }
    }
}

int main(){
    sc21::SC_input();

    sa(sc21::C);

    /*
    提出時にはここを消せ！！！！！！！！！！！
    std::cout<<simulator(sc21::C,sc21::I_PROB)<<"\n";
    */

    sc21::SC_output();
}