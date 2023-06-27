#ifndef _C_S_A_
#define _C_S_A_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

#ifndef uint
#define uint unsigned int
#endif

#ifndef TG
#define TG 0.1
#endif
#ifndef TA
#define TA 0.9
#endif

//#define MAXITER 200
//#define TGEN 0.1

#ifndef END
#define END 99
#endif

struct Optimizer
{
    int id;          //id
    double *curSol;     //Solução atual     [* Dimensões]
    double *sol;        //Solução provável  [* Dimensões]
    //double *temp;       //Auxiliar de troca [* Dimensões]
    struct drand48_data buffer;
    double curCost;     //Custo atual
    double cost;        //Custo da solução provável
    double prob;
    double result;

};

typedef struct
{
    char step;
    int iter;      //Interação
    int maxIter;
    
    int numOpt;    //Número de Otimizadores
    int dim;       //Número de Dimensões

    double tgen;    //Temperatura de Geração
    double tac;     //Temperatura de Aceitação
    double gamma;
    double maxCost; //Valor máximo de custo
    double tmp;
    double prob;
    double probVar;

    double *bestSol;    //Melhor Solução [* Dimensões]
    double bestCost;    //Melhor valor de custo

    struct Optimizer *opts; //Otimizadores
    double **solution;      //Soluções parciais a serem retornadas

} CSA;

double ** csa_init(CSA * csa, int _numOpt, int _dim, int _maxIter); //, double **_points, int _num);
double ** csa_exec(CSA * csa, double *costs);
void csa_reset(CSA *csa);
void csa_finalize(CSA * csa);


#endif