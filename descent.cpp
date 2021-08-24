#include <bits/stdc++.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "./main.cpp"
namespace sc21 {
#include "sc21.h"
}

void change_c(){
    long double eps=1e-12;
    double C_back_1[N_GROUP][N_GROUP], C_back_2[N_GROUP][N_GROUP];
    double C_double[N_GROUP][N_GROUP];

    backprop(sc21::C, sc21::I_PROB, C_back_1);

    //勾配正なら0に 0なら維持 負なら1に
    //epsは10^-12(小さすぎる？)
    for(int i=0;i<N_GROUP;i++){
        for(int j=i+1;j<N_GROUP;j++){
            if(C_back_1[i][j]+C_back_1[j][i]>eps){
                sc21::C[i][j]=0;
                sc21::C[j][i]=0;
            }
            else if(C_back_1[i][j]+C_back_1[j][i]<-1*eps){
                sc21::C[i][j]=1;
                sc21::C[j][i]=1;
            }
        }
    }
}

void descent(){
    int q=100;
    while(q--){
        change_c();
    }
}

int main(){
    sc21::SC_input();

    //初期化(全部1で埋める)
    for(int i=0;i<N_GROUP;i++){
        sc21::C[i][i]=0;
        for(int j=i+1;j<N_GROUP;j++){
            sc21::C[i][j]=1;
            sc21::C[j][i]=1;
        }
    }

    descent();

    sc21::SC_output();
}