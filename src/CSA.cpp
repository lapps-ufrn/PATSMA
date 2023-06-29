#include "CSA.hpp"

#include <iostream>  // cout, endl
#include <limits>    // numeric_limits

#define DBL_MAX std::numeric_limits<double>::max()

inline void CSA::copy(double *out, double *in) const {
  for (int i = 0; i < dim; i++) {
    out[i] = in[i];
  }
}

/**
 * Make round shift for values < -1 and > 1
 * value : in - Point
 * return : out - Point between -1 and 1
 */
auto CSA::rotate(double value) -> double {
  int i = (int)value;
  if (value > 1.0) {
    return (-1 + (value - i));
  } else if (value < -1.0) {
    return (1 + (value - i));
  }
  return value;
}

/**
 * Find the cost maximum value and save in maxCost
 */
void CSA::maxCost() {
  max_cost = opts[0].curCost;
  for (int k = 1; k < num_opt; k++) {
    if (opts[k].curCost > max_cost) {
      max_cost = opts[k].curCost;
    }
  }
}

/**
 * Switch values in vector position [i] from current solution to solution,
 * same from current cost to cost and check if this new cost is the maximum
 * i : in - Switch position
 */
void CSA::swap(int i) {
  double *temp;
  tmp = opts[i].cost;
  opts[i].cost = opts[i].curCost;
  opts[i].curCost = tmp;

  temp = opts[i].sol;
  opts[i].sol = opts[i].curSol;
  opts[i].curSol = temp;

  if (opts[i].curCost > max_cost) {
    max_cost = opts[i].curCost;
  }
}

/**
 * Reset the CSA
 * level : in - Select the level of reseting
 *    level 2 - Reset the number of iteretions
 *    level 1 - Reset the points and the temperatures (plus the previous ones)
 *    level 0 - Remove the best solution (plus the previous ones)
 *
 */
void CSA::reset(int level) {
  int i, j;
  switch (level) {
    case 0:  // Remove the infomation of best solution
      best_cost = DBL_MAX;

    case 1:  // Reset the points and the csa
      step = 0;
      tgen = TG;
      tac = TA;

      // Step 2.2: Inital points generate
      for (i = 0; i < num_opt; i++) {
        for (j = 0; j < dim; j++) {
          drand48_r(&opts[i].buffer, &opts[i].result);
          opts[i].curSol[j] = (opts[i].result * 2.0 - 1.0);
          solution[i][j] = opts[i].curSol[j];
        }
        opts[i].curCost = 0;
      }

    case 2:  // Reset only the number of iteretions
      iter = 0;
      break;

    default:
      std::cout << "There are not the CSA reset option level " << level
                << std::endl;
      break;
  }
}

/**
 * Variables inicialization
 * _num_opt : in - Amount of optimizer
 * _dim : in - Cost Function Dimension
 * _max_iter : in - Maximun iteration
 */
CSA::CSA(int _num_opt, int _dim, int _max_iter) {
  int i, j;

  this->iter = 0;
  this->step = 0;
  this->num_opt = _num_opt;
  this->dim = _dim;
  try {
    this->max_iter = (int)(_max_iter / (double)_num_opt);
  } catch (std::runtime_error &e) {
    std::cout << "Exception occurred" << std::endl << e.what();
  }

  this->tgen = TG;
  this->tac = TA;
  this->gamma = 0.0;
  this->best_cost = DBL_MAX;

  this->max_cost = 0;
  this->tmp = 0;
  this->prob = 0;
  this->prob_var = 0;

  try {
    opts = new Opt[num_opt];
    solution = new double *[num_opt];
    best_sol = new double[dim];
  } catch (const std::bad_alloc &e) {
    std::cout << "Memory Allocation"
              << " is failed: " << e.what() << std::endl;
  }

  // Step 1: Initialize variables [Optimizers]
  srand(time(nullptr));
  try {
    for (i = 0; i < num_opt; i++) {
      opts[i].id = i;
      opts[i].curSol = new double[dim];
      opts[i].sol = new double[dim];
      solution[i] = new double[dim];
      srand48_r(rand(), &opts[i].buffer);
    }
  } catch (const std::bad_alloc &e) {
    std::cout << "Memory Allocation"
              << " is failed: " << e.what() << std::endl;
  }

  // Step 2.1: Inital points copy
  //  for (i=0; i<_num_opt && i < _num; i++)
  //  {
  //      for(j=0;j<_dim;j++)
  //      {
  //          opts[i].curSol[j] = _points[i][j];
  //      }
  //      solution[i] = opts[i].curSol;
  //      copy(solution[i], opts[i].curSol, dim );
  //  }

  // Step 2.2: Inital points generate
  for (i = 0; i < num_opt; i++) {
    for (j = 0; j < dim; j++) {
      drand48_r(&opts[i].buffer, &opts[i].result);
      opts[i].curSol[j] = (opts[i].result * 2.0 - 1.0);
      solution[i][j] = opts[i].curSol[j];
    }
  }
}

/**
 * Coupled Simulated Annealing function
 * csa : in - Global variables
 * costs : in - Cost vector for all Optimizers
 * return : out - Points vector for all Optimizers
 */
void CSA::partial_exec(double *costs) {
  int k, j, i;

  while (step != END) {
    switch (step) {
      case 0:
        // Step 2.3: Find the best solution among the initial points and add
        // initial costs
        for (i = 0; i < num_opt; i++) {
          opts[i].curCost = costs[i];
          if (costs[i] < best_cost) {
            best_cost = costs[i];
            copy(best_sol, opts[i].curSol);
          }
        }
        // Step 3: Calculate gamma
        maxCost();
        gamma = 0;
        for (i = 0; i < num_opt; i++) {
          gamma += exp((opts[i].curCost - max_cost) / tac);
        }

      case 1:
        // Step 4: Generated New Points
        for (i = 0; i < num_opt; i++) {
          for (j = 0; j < dim; j++) {
            drand48_r(&opts[i].buffer, &opts[i].result);
            opts[i].result = tan(M_PI * (opts[i].result - 0.5));
            opts[i].sol[j] = rotate(opts[i].curSol[j] + tgen * opts[i].result);
          }
          copy(solution[i], opts[i].sol);
          // opts[i].sol[0] = solution[i][0];
          // for (int ii = 0; ii < dim; ii++) {
          //   opts[i].sol[ii] = solution[i][ii];
          // }
        }
        step = 2;
        return;

      case 2:
        // Step 5: Define if accept new solutions
        for (i = 0; i < num_opt; i++) {
          opts[i].cost = costs[i];

          // Step 5.1: If new soluiton is better
          if (opts[i].cost < opts[i].curCost) {
            // Step 5.1.2: Better global solution
            if (opts[i].cost < best_cost) {
              best_cost = opts[i].cost;
              copy(best_sol, opts[i].sol);
            }
            swap(i);
          }
          // Step 5.2: Else test probability of accept
          else {
            drand48_r(&opts[i].buffer, &opts[i].result);
            opts[i].prob = exp((opts[i].curCost - max_cost) / tac) / gamma;

            if (opts[i].prob > opts[i].result) {
              swap(i);
            }
          }
        }  // Optimizers

        // Step 6: Stop criterium
        //  if((++stab_iter) >= max_iter_stab)
        //  {
        //      step = END;
        //      break;
        //  }
        //  else
        if ((++iter) >= max_iter) {
          step = END;
          break;
        }
        // Step 7: Update variables
        else {
          step = 1;
          // Step 7.1: Procura de Máximo
          // maxCost(csa);//Esse teste já está sendo feito do swap()
          // Step 7.2: Calculate gamma (same that step 3)
          gamma = tmp = 0.0;
          for (k = 0; k < num_opt; k++) {
            gamma += exp((opts[k].curCost - max_cost) / tac);
            tmp += exp(2.0 * (opts[k].curCost - max_cost) / tac);
          }
          // Step 7.3: Update accept temperature
          tmp = tmp / (gamma * gamma);
          prob_var =
              (tmp * ((double)(num_opt)) - 1.0) / ((double)num_opt - 1.0);
          if (prob_var >= 0.99) {
            tac += 0.05 * tac;
          } else {
            tac -= 0.05 * tac;
          }
          // Step 7.5: Update generation temperature
          // tgen = tgen * log((double)iter+2.0)/log((double)iter+3.0);
          // tgen = tgen * ((double)iter+1.0)/((double)iter+2.0);
          // tgen = TGEN/log((double)iter+2.);
          // tgen = TGEN/((double)iter+1.);
          tgen = 0.99997 * tgen;
        }
        break;
    }
  }
}
