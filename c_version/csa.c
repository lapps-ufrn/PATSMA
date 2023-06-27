#include "csa.h"

#define DBL_MAX 1.7976931348623158e+308
#define copy(a, b, s) \
  for (int ii = 0; ii < s; ii++) b[ii] = a[ii];
#define reScaleToInt(in, out, min, max, size) \
  for (int ii = 0; ii < size; ii++)           \
    out[ii] = ((in[ii] + 1) / 2) * (max - min) + min;

#define debug printf

/**
 * Make round shift for values < -1 and > 1
 * value : in - Point
 * return : out - Point between -1 and 1
 */
double rotate(double value) {
  int i = value;
  if (value > 1.0)
    return (-1 + (value - i));
  else if (value < -1.0)
    return (1 + (value - i));
  return value;
}

/**
 * Find the cost maximum value and save in csa->maxCost
 * csa : in - Global variables
 */
void maxCost(CSA *csa) {
  csa->maxCost = csa->opts[0].curCost;
  for (int k = 1; k < csa->numOpt; k++) {
    if (csa->opts[k].curCost > csa->maxCost) {
      csa->maxCost = csa->opts[k].curCost;
    }
  }
}

/**
 * Switch values in vector position [i] from current solution to solution,
 * same from current cost to cost and check if this new cost is the maximum
 * csa : in - Global variables
 * i : in - Switch position
 */
void swap(CSA *csa, int i) {
  double *temp;
  csa->tmp = csa->opts[i].cost;
  csa->opts[i].cost = csa->opts[i].curCost;
  csa->opts[i].curCost = csa->tmp;

  temp = csa->opts[i].sol;
  csa->opts[i].sol = csa->opts[i].curSol;
  csa->opts[i].curSol = temp;

  if (csa->opts[i].curCost > csa->maxCost) {
    csa->maxCost = csa->opts[i].curCost;
  }
}

void csa_reset(CSA *csa) {
  csa->iter = 0;
  csa->step = 0;
  csa->tgen = TG;
  csa->tac = TA;
  csa->gamma = 0.0;

  for (i = 0; i < _numOpt; i++) {
    for (j = 0; j < _dim; j++) {
      drand48_r(&csa->opts[i].buffer, &csa->opts[i].result);
      csa->opts[i].curSol[j] = (csa->opts[i].result * 2.0 - 1.0);
    }
    csa->solution[i] = csa->opts[i].curSol;
    copy(csa->opts[i].curSol, csa->solution[i], csa->dim);
  }
}

/**
 * Variables inicialization
 * csa : in - Global variables
 * _numOpt : in - Amount of optimizer
 * _dim : in - Cost Function Dimension
 * _maxIter : in - Maximun iteration
 * return : out - Initial points
 */
double **csa_init(CSA *csa, int _numOpt, int _dim, int _maxIter) {
  csa->iter = 0;
  csa->step = 0;
  csa->numOpt = _numOpt;
  csa->dim = _dim;
  csa->maxIter = _maxIter / _numOpt;

  csa->tgen = TG;
  csa->tac = TA;
  csa->gamma = 0.0;
  csa->bestCost = DBL_MAX;

  csa->opts = (struct Optimizer *)malloc(_numOpt * sizeof(struct Optimizer));
  if (csa->opts == NULL) {
    printf("Failed allocating memory for Optimizer data\n");
    exit(0);
  }
  csa->solution = (double **)malloc(_numOpt * sizeof(double *));
  if (csa->solution == NULL) {
    printf("Failed allocating memory for global solution data\n");
    exit(0);
  }
  csa->bestSol = (double *)malloc(_dim * sizeof(double));
  if (csa->bestSol == NULL) {
    printf("Failed allocating memory for best solution data\n");
    exit(0);
    ;
  }

  int i, j;
  // Step 1: Initialize variables [Optimizers]
  srand(time(NULL));
  for (i = 0; i < _numOpt; i++) {
    csa->opts[i].id = i;
    csa->opts[i].curSol = (double *)malloc(_dim * sizeof(double));
    if (csa->opts[i].curSol == NULL) {
      printf(
          "Failed allocating memory for Optimizer Current solution data [%i]\n",
          i);
      exit(0);
    }
    csa->opts[i].sol = (double *)malloc(_dim * sizeof(double));
    if (csa->opts[i].sol == NULL) {
      printf("Failed allocating memory for Optimizer solution data [%i]\n", i);
      exit(0);
    }
    csa->solution[i] = (double *)malloc(_dim * sizeof(double));
    if (csa->solution[i] == NULL) {
      printf("Failed allocating memory for global solution data [%i]\n", i);
      exit(0);
    }
    srand48_r(rand(), &csa->opts[i].buffer);
  }

  // Step 2.1: Inital points copy
  //  for (i=0; i<_numOpt && i < _num; i++)
  //  {
  //      for(j=0;j<_dim;j++)
  //      {
  //          csa->opts[i].curSol[j] = _points[i][j];
  //      }
  //      csa->solution[i] = csa->opts[i].curSol;
  //      copy( csa->opts[i].curSol, csa->solution[i], csa->dim );
  //  }

  // Step 2.2: Inital points generate
  for (i = 0; i < _numOpt; i++) {
    for (j = 0; j < _dim; j++) {
      drand48_r(&csa->opts[i].buffer, &csa->opts[i].result);
      csa->opts[i].curSol[j] = (csa->opts[i].result * 2.0 - 1.0);
    }
    csa->solution[i] = csa->opts[i].curSol;
    copy(csa->opts[i].curSol, csa->solution[i], csa->dim);
  }

  return csa->solution;
}

/**
 * Coupled Simulated Annealing function
 * csa : in - Global variables
 * costs : in - Cost vector for all Optimizers
 * return : out - Points vector for all Optimizers
 */
double **csa_exec(CSA *csa, double *costs) {
  int k, j, i;

  while (csa->step != END) {
    switch (csa->step) {
      case 0:
        csa->bestCost = csa->opts[0].curCost = costs[0];
        copy(csa->opts[0].curSol, csa->bestSol, csa->dim);

        for (i = 1; i < csa->numOpt; i++) {
          csa->opts[i].curCost = costs[i];
          if (costs[i] < csa->bestCost) {
            csa->bestCost = costs[i];
            copy(csa->opts[i].curSol, csa->bestSol, csa->dim);
          }
        }
        // Step 3: Calculate gexit(0);amma
        maxCost(csa);
        for (i = 0; i < csa->numOpt; i++) {
          csa->gamma += exp((csa->opts[i].curCost - csa->maxCost) / csa->tac);
        }

      case 1:
        // Step 4: Generated New Points
        for (i = 0; i < csa->numOpt; i++) {
          for (j = 0; j < csa->dim; j++) {
            drand48_r(&csa->opts[i].buffer, &csa->opts[i].result);
            csa->opts[i].result = tan(PI * (csa->opts[i].result - 0.5));
            csa->opts[i].sol[j] = rotate(csa->opts[i].curSol[j] +
                                         csa->tgen * csa->opts[i].result);
          }
          copy(csa->opts[i].sol, csa->solution[i], csa->dim);
        }
        csa->step = 2;
        return csa->solution;

      case 2:
        // Step 5: Define if accept new solutions
        for (i = 0; i < csa->numOpt; i++) {
          csa->opts[i].cost = costs[i];

          // Step 5.1: If new soluiton is better
          if (csa->opts[i].cost < csa->opts[i].curCost) {
            // Step 5.1.2: Better global solution
            if (csa->opts[i].cost < csa->bestCost) {
              csa->bestCost = csa->opts[i].cost;
              copy(csa->opts[i].sol, csa->bestSol, csa->dim);
            }
            swap(csa, i);
          }
          // Step 5.2: Else test probability of accept
          else {
            drand48_r(&csa->opts[i].buffer, &csa->opts[i].result);
            csa->opts[i].prob =
                exp((csa->opts[i].curCost - csa->maxCost) / csa->tac) /
                csa->gamma;

            if (csa->opts[i].prob > csa->opts[i].result) {
              swap(csa, i);
            }
          }
        }  // Optimizers

        // Step 6: Stop criterium
        //  if((++csa->stab_iter) >= csa->maxIter_stab)
        //  {
        //      csa->step = END;
        //      break;
        //  }
        //  else
        if ((++csa->iter) >= csa->maxIter) {
          csa->step = END;
          break;
        }
        // Step 7: Update variables
        else {
          csa->step = 1;
          // Step 7.1: Procura de Máximo
          // maxCost(csa);//Esse teste já está sendo feito do swap()
          // Step 7.2: Calculate gamma (same that step 3)
          csa->gamma = csa->tmp = 0.0;
          for (k = 0; k < csa->numOpt; k++) {
            csa->gamma += exp((csa->opts[k].curCost - csa->maxCost) / csa->tac);
            csa->tmp +=
                exp(2.0 * (csa->opts[k].curCost - csa->maxCost) / csa->tac);
          }
          // Step 7.3: Update accept temperature
          csa->tmp = csa->tmp / (csa->gamma * csa->gamma);
          csa->probVar = (csa->tmp * ((double)(csa->numOpt)) - 1.0) /
                         ((double)csa->numOpt - 1.0);
          if (csa->probVar >= 0.99)
            csa->tac += 0.05 * csa->tac;
          else
            csa->tac -= 0.05 * csa->tac;
          // Step 7.5: Update generation temperature
          // csa->tgen = csa->tgen *
          // log((double)csa->iter+2.0)/log((double)csa->iter+3.0); csa->tgen =
          // csa->tgen * ((double)csa->iter+1.0)/((double)csa->iter+2.0);
          // csa->tgen = TGEN/log((double)csa->iter+2.);
          // csa->tgen = TGEN/((double)csa->iter+1.);
          csa->tgen = 0.99997 * csa->tgen;
        }
        break;
    }
  }
  return NULL;
}

/**
 * Memory free
 * csa : in - Global variables
 */
void csa_finalize(CSA *csa) {
  for (int i = 0; i < csa->numOpt; i++) {
    free(csa->opts[i].curSol);
    free(csa->opts[i].sol);
  }
  free(csa->opts);
  free(csa->solution);
}