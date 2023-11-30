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
  memcpy(out, in, m_dim * sizeof(double));
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
  double *temp = m_opts[i].probSol;
  m_opts[i].probSol = m_opts[i].curSol;
  m_opts[i].curSol = temp;

  m_tmp = m_opts[i].probCost;
  m_opts[i].probCost = m_opts[i].curCost;
  m_opts[i].curCost = m_tmp;

  if (m_opts[i].curCost > m_maxCost) {
    m_maxCost = m_opts[i].curCost;
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
      m_bestCost = DBL_MAX;
      memset(m_bestSol, 0, m_dim * sizeof(double));

    case 1:  // Reset points and temperatures
      m_step = 0;
      m_tGen = TG;
      m_tAcc = TA;
      m_maxCost = DBL_MIN;
      m_tmp = 0;
      m_probVar = 0;

      // Inital points generate
      for (i = 0; i < m_nOpt; i++) {
        for (j = 0; j < m_dim; j++) {
          drand48_r(&m_opts[i].buffer, &m_opts[i].result);
          m_opts[i].probSol[j] = (m_opts[i].result * 2.0 - 1.0);
          m_opts[i].curSol[j] = m_opts[i].probSol[j];  // Copy Initial Solutions
        }
        m_opts[i].curCost = DBL_MAX;
        m_opts[i].probCost = DBL_MAX;
      }

    case 2:  // Reset only number of iteretions
      m_iter = 0;
      m_iOpt = 0;
      if (m_step == END) m_step = 1;
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

  this->m_iOpt = 0;
  this->m_iter = 0;
  this->m_step = 0;
  this->m_nOpt = _num_opt;
  this->m_dim = _dim;
  try {
    this->m_maxIter = (int)(_max_iter / (double)_num_opt);
  } catch (std::runtime_error &e) {
    std::cout << "Exception occurred" << std::endl << e.what();
  }

  this->m_tGen = TG;
  this->m_tAcc = TA;
  this->m_gamma = 0.0;
  this->m_bestCost = DBL_MAX;

  this->m_maxCost = DBL_MIN;
  this->m_tmp = 0;
  this->m_probVar = 0;

  try {
    m_point = new double[m_dim];
    m_opts = new Opt[m_nOpt];
    m_bestSol = new double[m_dim];
    for (i = 0; i < m_nOpt; i++) {
      m_opts[i].curSol = new double[m_dim];
      m_opts[i].probSol = new double[m_dim];
    }
  } catch (const std::bad_alloc &e) {
    std::cout << "Memory Allocation"
              << " is failed: " << e.what() << std::endl;
  }

  // Step 1: Initialize variables [Optimizers]
  srand(time(nullptr));
  for (i = 0; i < m_nOpt; i++) {
    m_opts[i].id = i;
    m_opts[i].curCost = DBL_MAX;
    m_opts[i].probCost = DBL_MAX;
    srand48_r(rand(), &m_opts[i].buffer);
  }

  // Step 1.1: Generate inital points
  for (i = 0; i < m_nOpt; i++) {
    for (j = 0; j < m_dim; j++) {
      drand48_r(&m_opts[i].buffer, &m_opts[i].result);
      m_opts[i].probSol[j] = (m_opts[i].result * 2.0 - 1.0);
      m_opts[i].curSol[j] = m_opts[i].probSol[j];  // Copy Initial Solutions
    }
  }
}

/// @brief Coupled Simulated Annealing function
/// @param costs Cost vector for all Optimizers
void CSA::partial_exec() {
  int k, j, i;

  while (m_step != END) {
    switch (m_step) {
      case 0:
        // Step 2: Initial points
        for (i = 0; i < m_nOpt; i++) {
          m_opts[i].curCost = m_opts[i].probCost;  // Copy initial costs
          if (m_opts[i].curCost < m_bestCost) {    // Find the best solution
            m_bestCost = m_opts[i].curCost;
            copy_solution(m_bestSol, m_opts[i].curSol);
          }
          if (m_maxCost < m_opts[i].curCost) {  // Find the max cost
            m_maxCost = m_opts[i].curCost;
          }
        }
        // Step 3: Calculate m_gamma
        m_gamma = 0;
        for (i = 0; i < m_nOpt; i++) {
          m_gamma += exp((m_opts[i].curCost - m_maxCost) / m_tAcc);
        }

      case 1:
        // Step 4: Generated New Points
        for (i = 0; i < m_nOpt; i++) {
          for (j = 0; j < m_dim; j++) {
            drand48_r(&m_opts[i].buffer, &m_opts[i].result);
            m_opts[i].result = tan(M_PI * (m_opts[i].result - 0.5));
            m_opts[i].probSol[j] =
                rotate(m_opts[i].curSol[j] + m_tGen * m_opts[i].result);
          }
        }
        m_step = 2;
        return;

      case 2:
        // Step 5: Define if accept new solution
        for (i = 0; i < m_nOpt; i++) {
          // Step 5.1: If new soluiton is better
          if (m_opts[i].probCost < m_opts[i].curCost) {
            swap_opt_info(i);
            // Step 5.1.2: Better global solution
            if (m_opts[i].curCost < m_bestCost) {
              m_bestCost = m_opts[i].curCost;
              copy_solution(m_bestSol, m_opts[i].curSol);
            }
          }
          // Step 5.2: Else test probability of accept
          else {
            drand48_r(&m_opts[i].buffer, &m_opts[i].result);
            m_opts[i].prob =
                exp((m_opts[i].curCost - m_maxCost) / m_tAcc) / m_gamma;

            if (m_opts[i].prob > m_opts[i].result) {
              swap_opt_info(i);
            }
          }
        }  // Optimizers

        // Step 6: Stop criterium
        if ((++m_iter) >= m_maxIter) {
          m_step = END;
          break;
        }
        // Step 7: Update variables
        else {
          m_step = 1;
          // Step 7.1: Calculate m_gamma (same that m_step 3)
          m_gamma = m_tmp = 0.0;
          for (k = 0; k < m_nOpt; k++) {
            m_gamma += exp((m_opts[k].curCost - m_maxCost) / m_tAcc);
            m_tmp += exp(2.0 * (m_opts[k].curCost - m_maxCost) / m_tAcc);
          }
          // Step 7.2: Update accept temperature
          m_tmp = m_tmp / (m_gamma * m_gamma);
          m_probVar =
              (m_tmp * ((double)(m_nOpt)) - 1.0) / ((double)m_nOpt - 1.0);
          if (m_probVar >= 0.99) {
            m_tAcc += 0.05 * m_tAcc;
          } else {
            m_tAcc -= 0.05 * m_tAcc;
          }
          // Step 7.3: Update generation temperature
          // m_tGen = m_tGen * log((double)m_iter+2.0)/log((double)m_iter+3.0);
          // m_tGen = m_tGen * ((double)m_iter+1.0)/((double)m_iter+2.0);
          // m_tGen = TGEN/log((double)m_iter+2.);
          // m_tGen = TGEN/((double)m_iter+1.);
          m_tGen = 0.99997 * m_tGen;
        }
        break;
    }
  }
}

double *CSA::run(double cost) {
  // Save costs as Probable costs values
  if (m_iOpt > 0 && m_iOpt <= m_nOpt) {
    m_opts[m_iOpt - 1].probCost = cost;
  }
  // Return to aplication the next solution
  if (m_iOpt < m_nOpt) {
    copy_solution(m_point, m_opts[m_iOpt].probSol);
    m_iOpt += 1;
  } else {
    partial_exec();
    if (m_step == END) {  // End execution, return best solution
      copy_solution(m_point, m_bestSol);
    } else {  // Return first optimazate
      copy_solution(m_point, m_opts[0].probSol);
      m_iOpt = 1;
    }
  }

  return m_point;
}

/// @brief Destructor
CSA::~CSA() {
  delete[] m_point;
  delete[] m_bestSol;
  for (int i = 0; i < m_nOpt; i++) {
    delete[] m_opts[i].curSol;
    delete[] m_opts[i].probSol;
  }
  delete[] m_opts;
}

void CSA::print() const {
  std::cout << "----------- CSA Parameters -----------" << std::endl;
  std::cout << "Dim: " << m_dim << "\tNOpt: " << m_nOpt
            << "\tNEval: " << (m_maxIter * m_nOpt) << std::endl;
  std::cout << "{";
  for (int i = 0; i < m_nOpt; i++) {
    std::cout << "[" << i << ", ";
    std::cout << "(";
    for (int j = 0; j < m_dim; j++) {
      std::cout << m_opts[i].probSol[j] << ((j < m_dim - 1) ? "," : "");
    }
    std::cout << ")";
    std::cout << "]" << ((i < m_nOpt - 1) ? "\t" : "");
  }
  std::cout << "}" << std::endl;
  std::cout << "---------------------------------------" << std::endl;
}