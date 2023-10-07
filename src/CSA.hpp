#ifndef _CSA_
#define _CSA_

#include <cmath>
#include <ctime>  // drand48_data

#include "NumericalOptimizer.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
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

// #define MAXITER 200
// #define TGEN 0.1

#define END 0x99

class CSA : public NumericalOptimizer {
  struct Opt {
    int id;           // id
    double *curSol;   // Solução atual     [* Dimensões]
    double *probSol;  // Solução provável  [* Dimensões]
    struct drand48_data buffer;
    double curCost;   // Custo atual
    double probCost;  // Custo da solução provável
    // Auxiliars
    double prob;
    double result;
  };

  int step;
  int iter;      // Iteration
  int max_iter;  // Max number of iterations

  int i_opt;    // Iterator for Optimizers
  int num_opt;  // Number of Optimizers
  int dim;      // Number of Dimensions

  double tgen;  // Generation Temperature
  double tac;   // Acceptance Temperature
  double gamma;
  double max_cost;  // Maximum cost value
  double tmp;
  double prob_var;

  double *best_sol;  // Best Solution [* Dimensões]
  double best_cost;  // Best Cost relative to best solution

  struct Opt *opts;  // Optimizers
  double *point;     // Point to return

  void swap_opt_info(int i);
  void copy_solution(double *out, double *in) const;
  static double rotate(double value);
  void partial_exec();

  CSA() = delete;
  CSA operator=(CSA) = delete;
  CSA &operator=(CSA &&) = delete;
  CSA(const CSA &) = delete;
  CSA(CSA &&) = delete;

 public:
  int getNumPoints() const override { return num_opt; }
  int getDimension() const override { return dim; }
  void print() const override;
  void reset(int level) override;
  bool isEnd() const override { return step == END; }
  double *run(double cost) override;

  CSA(int _num_opt, int _dim, int _max_iter);
  ~CSA();
};

#endif