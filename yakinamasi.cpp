#include <bits/stdc++.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "./main.cpp"

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

// 焼きなまし法
void sa(int C[][N_GROUP]){
    double start_time=sc21::TIME0;

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
    double TIME_LIMIT=295.0;

    while(true){
        if(sc21::TIME0-start_time>TIME_LIMIT)break;
        auto new_state=C;

        //だんだん変更量を減少させるように
        int change=1+10-10*(sc21::TIME0-start_time)/TIME_LIMIT;
        
        modify(new_state,change);
        double new_score=simulator(new_state,sc21::I_PROB);
        double pre_score=simulator(C,sc21::I_PROB);

        double temp=start_temp+(end_temp-start_temp)*(sc21::TIME0-start_time)/TIME_LIMIT;
        double prob=exp((new_score-pre_score)/temp);

        if(prob>(xor128()%1280000)/1280000){
            C=new_state;
        }
    }
}