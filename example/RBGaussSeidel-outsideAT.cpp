#include <omp.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "Autotuning.hpp"

#define N 2001
#define MAX_ITER 1000

int PARAM_AMOUNT;

void initialize(double **A, int n) {
  int i, j;
  for (j = 0; j < n + 1; j++) {
    A[0][j] = 1.0;
  }
  for (i = 1; i < n + 1; i++) {
    A[i][0] = 1.0;
    for (j = 1; j < n + 1; j++) {
      A[i][j] = 0.0;
    }
  }
}

double matrix_calculation(double **A, int n, int *chunk) {
  double tmp;
  double diff = 0;
  int i, j;

#ifdef _OPENMP
#pragma omp parallel private(tmp, i, j)
  {
#pragma omp for reduction(+ : diff) schedule(dynamic, chunk[0])
#endif
    for (i = 1; i <= n; ++i) {
      for (j = 1; j <= n; ++j) {
        if ((i + j) % 2 == 1) {
          tmp = A[i][j];
          A[i][j] = 0.2 * (A[i][j] + A[i][j - 1] + A[i - 1][j] + A[i][j + 1] + A[i + 1][j]);
          diff += fabs(A[i][j] - tmp);
        }
      }
    }
#ifdef _OPENMP
#pragma omp for reduction(+ : diff) schedule(dynamic, (PARAM_AMOUNT == 1 ? chunk[0] : chunk[1]))
#endif
    for (i = 1; i <= n; ++i) {
      for (j = 1; j <= n; ++j) {
        if ((i + j) % 2 == 0) {
          tmp = A[i][j];
          A[i][j] = 0.2 * (A[i][j] + A[i][j - 1] + A[i - 1][j] + A[i][j + 1] + A[i + 1][j]);
          diff += fabs(A[i][j] - tmp);
        }
      }
    }
#ifdef _OPENMP
  }
#endif
  return diff;
}

void solve_parallel(double **A, int n) {
  /*PATSMA parameters*/
  int dim = PARAM_AMOUNT;
  const int min = 1;
  const int max = (N - 2) / (omp_get_max_threads() * 2);
  const int ignore = 0;
  const int n_opt = 4;
  const int n_iter = 20;

  /*PATSMA INSTANCE*/
  int *chunk = new int[dim];
  auto *at = new Autotuning(min, max, ignore, dim, n_opt, n_iter);
  int iters;
  double diff = 0;

  printf("\n\n-----------------------Parallel Red Black Solver-----------------------\n\n\n");

  at->entireExecRuntime(matrix_calculation, chunk, A, N - 1);
  initialize(A, N);

  for (iters = 1; iters < MAX_ITER; ++iters) {
    diff = matrix_calculation(A, N - 1, chunk);
  }

  printf("Difference after %3d iterations: %f\n", iters, diff);
  printf("\n\nIteration LIMIT Reached...Exiting\n\n");

  if (PARAM_AMOUNT == 1) {
    printf("PATSMA calculed point: %d\n", chunk[0]);
  } else {
    printf("PATSMA calculed points: %d %d\n", chunk[0], chunk[1]);
  }
  at->print();
  delete at;
}

int main(int argc, char *argv[]) {

  PARAM_AMOUNT = 1;
  if (argc > 1) {
    PARAM_AMOUNT = atoi(argv[1]);
  }

  int i;
  double t_start, t_end;
  double **A;
  A = new double *[N + 2];
  for (i = 0; i < N + 2; i++) {
    A[i] = new double[N + 2];
  }

  initialize(A, N);

  t_start = omp_get_wtime();
  solve_parallel(A, N - 1);
  t_end = omp_get_wtime();

  printf("Computation time for Parallel approach(secs) on %i threads: %f\n\n",
         omp_get_max_threads(), t_end - t_start);
}