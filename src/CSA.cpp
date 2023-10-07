#include "CSA.hpp"

#include <cstring>   // memset, memcpy
#include <iostream>  // cout, endl

#ifndef DBL_MAX
#include <limits>  // numeric_limits
#define DBL_MAX std::numeric_limits<double>::max()
#define DBL_MIN std::numeric_limits<double>::min()
#endif

/// @brief Copy solution vector
/// @param out the output solution vector
/// @param in the input solution vector
inline void CSA::copy_solution(double *out, double *in) const {
  memcpy(out, in, dim * sizeof(double));
}

/// @brief Make round shift for values < -1 and > 1
/// @param value Point
/// @return Point between -1 and 1
auto CSA::rotate(double value) -> double {
  int i = (int)value;
  if (value > 1.0) {
    return (-1 + (value - i));
  } else if (value < -1.0) {
    return (1 + (value - i));
  }
  return value;
}

/// @brief Switch values in vector position [i] from current solution to
/// solution, same from current cost to cost and check if this new cost is the
/// maximum
/// @param i Switch position
void CSA::swap_opt_info(int i) {
  double *temp = opts[i].probSol;
  opts[i].probSol = opts[i].curSol;
  opts[i].curSol = temp;

  tmp = opts[i].probCost;
  opts[i].probCost = opts[i].curCost;
  opts[i].curCost = tmp;

  if (opts[i].curCost > max_cost) {
    max_cost = opts[i].curCost;
  }
}

/// @brief Reset the CSA
/// @param level Select the level of reseting
///    level 2 - Reset the number of iteretions
///    level 1 - Reset the points and the temperatures (plus the previous ones)
///    level 0 - Remove the best solution (plus the previous ones)
void CSA::reset(int level) {
  int i, j;
  switch (level) {
    case 0:  // Remove infomation of best solution
      best_cost = DBL_MAX;
      memset(best_sol, 0, dim * sizeof(double));

    case 1:  // Reset points and temperatures
      step = 0;
      tgen = TG;
      tac = TA;
      max_cost = DBL_MIN;
      tmp = 0;
      prob_var = 0;

      // Inital points generate
      for (i = 0; i < num_opt; i++) {
        for (j = 0; j < dim; j++) {
          drand48_r(&opts[i].buffer, &opts[i].result);
          opts[i].probSol[j] = (opts[i].result * 2.0 - 1.0);
          opts[i].curSol[j] = opts[i].probSol[j];  // Copy Initial Solutions
        }
        opts[i].curCost = DBL_MAX;
        opts[i].probCost = DBL_MAX;
      }

    case 2:  // Reset only number of iteretions
      iter = 0;
      i_opt = 0;
      if (step == END) step = 1;
      break;

    default:
      throw std::runtime_error("There is not the CSA reset option level " +
                               std::to_string(level));
      break;
  }
}

/// @brief Variables inicialization
/// @param _num_opt Amount of optimizer
/// @param _dim Cost Function Dimension
/// @param _max_iter Maximun iteration
CSA::CSA(int _num_opt, int _dim, int _max_iter) {
  if (_dim < 1) {
    throw std::invalid_argument("Dimensional Value Invalid! Set _dim > 0.");
  }
  if (_num_opt < 1) {
    throw std::invalid_argument("Optmizers Number Invalid! Set _num_opt > 0.");
  }
  if (_max_iter < 1) {
    throw std::invalid_argument(
        "Max number of intereration Invalid! Set _max_iter > 0.");
  }

  int i, j;

  this->i_opt = 0;
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

  this->max_cost = DBL_MIN;
  this->tmp = 0;
  this->prob_var = 0;

  try {
    point = new double[dim];
    opts = new Opt[num_opt];
    best_sol = new double[dim];
    for (i = 0; i < num_opt; i++) {
      opts[i].curSol = new double[dim];
      opts[i].probSol = new double[dim];
    }
  } catch (const std::bad_alloc &e) {
    std::cout << "Memory Allocation"
              << " is failed: " << e.what() << std::endl;
  }

  // Step 1: Initialize variables [Optimizers]
  srand(time(nullptr));
  for (i = 0; i < num_opt; i++) {
    opts[i].id = i;
    opts[i].curCost = DBL_MAX;
    opts[i].probCost = DBL_MAX;
    srand48_r(rand(), &opts[i].buffer);
  }

  // Step 1.1: Generate inital points
  for (i = 0; i < num_opt; i++) {
    for (j = 0; j < dim; j++) {
      drand48_r(&opts[i].buffer, &opts[i].result);
      opts[i].probSol[j] = (opts[i].result * 2.0 - 1.0);
      opts[i].curSol[j] = opts[i].probSol[j];  // Copy Initial Solutions
    }
  }
}

/// @brief Coupled Simulated Annealing function
/// @param costs Cost vector for all Optimizers
void CSA::partial_exec() {
  int k, j, i;

  while (step != END) {
    switch (step) {
      case 0:
        // Step 2: Initial points
        for (i = 0; i < num_opt; i++) {
          opts[i].curCost = opts[i].probCost;  // Copy initial costs
          if (opts[i].curCost < best_cost) {   // Find the best solution
            best_cost = opts[i].curCost;
            copy_solution(best_sol, opts[i].curSol);
          }
          if (max_cost < opts[i].curCost) {  // Find the max cost
            max_cost = opts[i].curCost;
          }
        }
        // Step 3: Calculate gamma
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
            opts[i].probSol[j] =
                rotate(opts[i].curSol[j] + tgen * opts[i].result);
          }
        }
        step = 2;
        return;

      case 2:
        // Step 5: Define if accept new solution
        for (i = 0; i < num_opt; i++) {
          // Step 5.1: If new soluiton is better
          if (opts[i].probCost < opts[i].curCost) {
            swap_opt_info(i);
            // Step 5.1.2: Better global solution
            if (opts[i].curCost < best_cost) {
              best_cost = opts[i].curCost;
              copy_solution(best_sol, opts[i].curSol);
            }
          }
          // Step 5.2: Else test probability of accept
          else {
            drand48_r(&opts[i].buffer, &opts[i].result);
            opts[i].prob = exp((opts[i].curCost - max_cost) / tac) / gamma;

            if (opts[i].prob > opts[i].result) {
              swap_opt_info(i);
            }
          }
        }  // Optimizers

        // Step 6: Stop criterium
        if ((++iter) >= max_iter) {
          step = END;
          break;
        }
        // Step 7: Update variables
        else {
          step = 1;
          // Step 7.1: Calculate gamma (same that step 3)
          gamma = tmp = 0.0;
          for (k = 0; k < num_opt; k++) {
            gamma += exp((opts[k].curCost - max_cost) / tac);
            tmp += exp(2.0 * (opts[k].curCost - max_cost) / tac);
          }
          // Step 7.2: Update accept temperature
          tmp = tmp / (gamma * gamma);
          prob_var =
              (tmp * ((double)(num_opt)) - 1.0) / ((double)num_opt - 1.0);
          if (prob_var >= 0.99) {
            tac += 0.05 * tac;
          } else {
            tac -= 0.05 * tac;
          }
          // Step 7.3: Update generation temperature
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

double *CSA::run(double cost) {
  // Save costs as Probable costs values
  if (i_opt > 0 && i_opt <= num_opt) {
    opts[i_opt - 1].probCost = cost;
  }
  // Return to aplication the next solution
  if (i_opt < num_opt) {
    copy_solution(point, opts[i_opt].probSol);
    i_opt += 1;
  } else {
    partial_exec();
    if (step == END) {  // End execution, return best solution
      copy_solution(point, best_sol);
    } else {  // Return first optimazate
      copy_solution(point, opts[0].probSol);
      i_opt = 1;
    }
  }

  return point;
}

/// @brief Destructor
CSA::~CSA() {
  delete[] point;
  delete[] best_sol;
  for (int i = 0; i < num_opt; i++) {
    delete[] opts[i].curSol;
    delete[] opts[i].probSol;
  }
  delete[] opts;
}

void CSA::print() const {
  std::cout << "----------- CSA Parameters -----------" << std::endl;
  std::cout << "Dim: " << dim << "\tNOpt: " << num_opt
            << "\tNEval: " << (max_iter * num_opt) << std::endl;
  std::cout << "{";
  for (int i = 0; i < num_opt; i++) {
    std::cout << "[" << i << ", ";
    std::cout << "(";
    for (int j = 0; j < dim; j++) {
      std::cout << opts[i].probSol[j] << ((j < dim - 1) ? "," : "");
    }
    std::cout << ")";
    std::cout << "]" << ((i < num_opt - 1) ? "\t" : "");
  }
  std::cout << "}" << std::endl;
  std::cout << "---------------------------------------" << std::endl;
}