#ifndef _CSA_
#define _CSA_

#include <cmath>
#include <ctime>

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

#ifndef END
#define END 99
#endif

struct Opt {
  int id;          // id
  double *curSol;  // Solução atual     [* Dimensões]
  double *sol;     // Solução provável  [* Dimensões]
  // double *temp;       //Auxiliar de troca [* Dimensões]
  struct drand48_data buffer;
  double curCost;  // Custo atual
  double cost;     // Custo da solução provável
  // Auxiliars
  double prob;
  double result;
};

class CSA {
 private:
  friend class Autotuning;

  int step;
  int iter;  // Interação
  int max_iter;

  int num_opt;  // Número de Otimizadores
  int dim;      // Número de Dimensões

  double tgen;  // Temperatura de Geração
  double tac;   // Temperatura de Aceitação
  double gamma;
  double max_cost;  // Valor máximo de custo
  double tmp;
  double prob;
  double prob_var;

  double *best_sol;  // Melhor Solução [* Dimensões]
  double best_cost;  // Melhor valor de custo

  struct Opt *opts;   // Otimizadores
  double **solution;  // Soluções parciais a serem retornadas

  void maxCost();
  void swap(int i);
  void copy(double *out, double *in) const;
  static auto rotate(double value) -> double;

 public:
  void partial_exec(double *costs);
  void reset(int level);
  auto getSolution(int i) const -> double * { return solution[i]; }

  auto operator=(CSA) -> CSA = delete;
  auto operator=(CSA &&) -> CSA & = delete;
  CSA(const CSA &) = delete;
  CSA(CSA &&) = delete;

  CSA() = delete;
  CSA(int _num_opt, int _dim, int _max_iter);
  ~CSA() {
    delete[] best_sol;
    for (int i = 0; i < num_opt; i++) {
      delete[] opts[i].curSol;
      delete[] opts[i].sol;
      delete[] solution[i];
    }
    delete[] opts;
    delete[] solution;
  }
};

#endif