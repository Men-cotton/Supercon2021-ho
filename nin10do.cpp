#include <stdio.h>
#include <math.h>
#include "cssl.h"
#include <stdlib.h>
#include <omp.h>
#include "mpi.h"
#include "sc21.h"

const int MAX_STEP_GREEDY = 10000;
const int MAX_SEC_GREEDY = 290;
const int MAX_PROCESS_COUNT = 48;
const int INIT_SEED = 12345;

void initialize_network(int N_LINK, unsigned rnd, int id);

double S[N_GROUP];
double I[N_GROUP];
double R[N_GROUP];

double S0[N_GROUP];
double I0[N_GROUP];
double R0[N_GROUP];

inline void simulation(const int dayCount);

inline bool swap_C_targets(const int i1, const int j1, const int i2, const int j2);

inline double calc_score(const double p1[N_GROUP], const double p2[N_GROUP]);

inline void update_rnd(unsigned &rnd)
{
    rnd = rnd ^ (rnd << 13);
    rnd = rnd ^ (rnd >> 17);
    rnd = rnd ^ (rnd << 5);
}

int main(int argc, char **argv)
{
    unsigned seed[MAX_PROCESS_COUNT];
    /*-- for MPI --*/
    MPI_Status status;

    MPI_Init(&argc, &argv);
    double start = omp_get_wtime();
    SC_input();

    int myId, processCount;
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    /* generate seed */
    if (myId == 0)
    {
        int icon;
        double dwork[8];
        double r[MAX_PROCESS_COUNT];

        int rnd = INIT_SEED;
        c_dm_vranu5(&rnd, r, processCount, (long)0, dwork, &icon);
        for (int i = 0; i < processCount; i++)
        {
            seed[i] = (unsigned)(r[i] * 2147483647);
            // printf("#seed[%d]=%d\n", i, seed[i]);
        };
    };

    MPI_Bcast(&seed, processCount, MPI_INT, 0, MPI_COMM_WORLD);

    /* inference simulation */

    initialize_network(N_LINK, seed[myId], myId);

    simulation(T);

    double bestScore = calc_score(I_PROB, I);

    unsigned rnd = seed[myId];
    double *scoreArray;
    scoreArray = (double *)malloc(processCount * sizeof(double));

    for (int steps = 0; true /*steps < MAX_STEP_GREEDY*/; steps++)
    {
        update_rnd(rnd);
        bool swaped = false;
        int i1, i2, j1, j2;

        while (!swaped)
        {
            update_rnd(rnd);
            i1 = rnd % N_GROUP, j1 = (rnd / N_GROUP) % N_GROUP;
            update_rnd(rnd);
            i2 = rnd % N_GROUP, j2 = (rnd / N_GROUP) % N_GROUP;

            swaped = swap_C_targets(i1, j1, i2, j2);
        };

        simulation(T);

        double currentScore = calc_score(I_PROB, I);
        double deltaScore = currentScore - bestScore;

        if (deltaScore > 0.0)
        {
            update_rnd(rnd);
            if ((rnd % 100 - 20 > steps / 2000 && deltaScore < 200 - (steps / 1000)) || rnd % 10000 > 9990 + steps / 18000)
                bestScore = currentScore;
            else
                swap_C_targets(i1, j1, i2, j2);
        }
        else
            bestScore = currentScore;

        if (steps % 1000 == 0)
        {
            double now = omp_get_wtime();
            double sec = now - start;
            if (sec > MAX_SEC_GREEDY)
                break;

            /*
            if (myId == 10)
                printf("%d %d %u %lf \n", myId, steps, rnd, bestScore);
            */
        }

        /*
        if (steps % 10000 == 999)
        {
            printf("#%d steps:%d\n", myId, steps);

            MPI_Gather(&bestScore, 1, MPI_DOUBLE, scoreArray, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (myId == 0)
            {
                double bestScore = scoreArray[0];
                for (int i = 0; i < processCount; i++)
                {
                    // printf("#score[%d]=%lf\n", i, scoreArray[i]);
                    if (bestScore > scoreArray[i])
                    {
                        bestScore = scoreArray[i];
                        bestIndex = i;
                    };
                };

                // printf("#bestScore,bestIndex= %lf %d\n", bestScore, bestIndex);
            };
        }
        */
    };

    /*
    int edgeCount = 0;
    for (int i = 0; i < N_GROUP; i++)
    {
        for (int j = i; j < N_GROUP; j++)
        {
            if (C[i][j] == 1)
                edgeCount++;
        };
    };
    printf("#myId,edgeCount,N_LINK=%d %d %d \n", myId, edgeCount, N_LINK);

    for (int i = 0; i < N_GROUP; i++)
        printf("%d %d %f %f\n", myId, i, I[i], I_prob[i]);
    */

    int bestIndex = 0;
    MPI_Gather(&bestScore, 1, MPI_DOUBLE, scoreArray, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myId == 0)
    {
        double bestScore = scoreArray[0];
        for (int i = 0; i < processCount; i++)
        {
            // printf("#score[%d]=%lf\n", i, scoreArray[i]);
            if (bestScore > scoreArray[i])
            {
                bestScore = scoreArray[i];
                bestIndex = i;
            };
        };

        // printf("#bestScore,bestIndex= %lf %d\n", bestScore, bestIndex);
    };
    MPI_Bcast(&bestIndex, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (bestIndex != 0)
    {
        int itag = 0;
        if (myId == bestIndex)
        {
            MPI_Send(C, (N_GROUP * N_GROUP), MPI_INT, 0, itag, MPI_COMM_WORLD);
        }
        else if (myId == 0)
        {
            MPI_Recv(C, (N_GROUP * N_GROUP), MPI_INT, bestIndex, itag, MPI_COMM_WORLD, &status);
        };
    };

    if (myId == 0)
    {
        SC_output();
        printf("Score: %f\n", scoreArray[bestIndex]);
    }

    MPI_Finalize();
    return 0;
};

void initialize_network(int N_LINK, unsigned rnd, int id)
{
    for (int i = 0; i < N_GROUP; i++)
    {
        for (int j = 0; j < N_GROUP; j++)
        {
            C[i][j] = 0;
        };
    };

    int edgeCount = 0;
    while (edgeCount < N_LINK)
    {
        update_rnd(rnd);
        int i = rnd % N_GROUP, j = (rnd / N_GROUP) % N_GROUP;
        if (C[i][j] == 0 && i != j)
        {
            edgeCount++;
            C[i][j] = 1, C[j][i] = 1;
        };
    };
};

inline void simulation(const int dayCount)
{
#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < N_GROUP; i++)
        {
            S0[i] = (double)N[i];
            I0[i] = 0.0;
            R0[i] = 0.0;
        };

        I0[0]++, S0[0]--;

        for (int istep = 0; istep < dayCount; istep++)
        {
#pragma omp for
            for (int i = 0; i < N_GROUP; i++)
            {

                double contact = BETA * I0[i];
                for (int j = 0; j < N_GROUP; j++)
                {
                    contact += BETA2 * C[i][j] * I0[j];
                };

                S[i] = S0[i] - contact * S0[i];
                I[i] = I0[i] + contact * S0[i] - GAMMA * I0[i];
                R[i] = R0[i] + GAMMA * I0[i];
            };

#pragma omp for
            for (int i = 0; i < N_GROUP; i++)
            {
                S0[i] = S[i];
                I0[i] = I[i];
                R0[i] = R[i];
            };
        };
    };
};

inline double calc_score(const double p1[N_GROUP], const double p2[N_GROUP])
{
    double score = 0.0;
#pragma omp parallel for reduction(+ \
								   : score)
    for (int i = 0; i < N_GROUP; i++)
    {
        score += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    };
    return score;
};

inline bool swap_C_targets(const int i1, const int j1, const int i2, const int j2)
{
    if (i1 != j1 && i2 != j2 && C[i1][j1] != C[i2][j2])
    {
        std::swap(C[i1][j1], C[i2][j2]);
        std::swap(C[j1][i1], C[j2][i2]);
        return true;
    }
    else
    {
        return false;
    }
};