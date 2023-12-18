#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "../src/Autotuning.hpp"

#define tolerance 0.000012
#define max_iterations 100

#define N 8001
#define THREAD_COUNT 6

double **A, diff;

void display(double **V, int n) {
  for (int i = 0; i < n + 1; ++i) {
    for (int k = 0; k < n + 1; ++k) {
      printf("%f ", V[i][k]);
    }
    printf("\n");
  }
  printf("\n");
}

void initialize(double **A, int n) {
  int i, j;
  for (j = 0; j < n + 1; j++) {
    A[0][j] = 1.0;
  }
  for (i = 1; i < n + 1; i++) {
    A[i][0] = 1.0;
    for (j = 1; j < n + 1; j++) A[i][j] = 0.0;
  }
}

long usecs(void) {
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec * 1000000 + t.tv_usec;
}

double matrix_calculation(double **A, int n, int *chunk) {
  double tmp;
  double diff;
  int i, j;

#ifdef TWO_PARAM
  int chk_01 = chunk[0];
  int chk_02 = chunk[1];
#else  // ONE_PARAM
  int chk_01 = chunk[0];
  int chk_02 = chunk[0];
#endif

#pragma omp parallel num_threads(THREAD_COUNT) private(tmp, i, j)
  {
#pragma omp for reduction(+ : diff) schedule(dynamic, chk_01)
    for (i = 1; i <= n; ++i) {
      for (j = 1; j <= n; ++j) {
        if ((i + j) % 2 == 1) {
          tmp = A[i][j];
          A[i][j] = 0.2 * (A[i][j] + A[i][j - 1] + A[i - 1][j] + A[i][j + 1] + A[i + 1][j]);
          diff += fabs(A[i][j] - tmp);
        }
      }
    }
#pragma omp for reduction(+ : diff) schedule(dynamic, chk_02)
    for (i = 1; i <= n; ++i) {
      for (j = 1; j <= n; ++j) {
        if ((i + j) % 2 == 0) {
          tmp = A[i][j];
          A[i][j] = 0.2 * (A[i][j] + A[i][j - 1] + A[i - 1][j] + A[i][j + 1] + A[i + 1][j]);
          diff += fabs(A[i][j] - tmp);
        }
      }
    }
  }
  return diff;
}

void solve_parallel(double **A, int n) {
  /*PATSMA parameters*/
#ifdef TWO_PARAM
  const int dim = 2;
#else
  const int dim = 1;
#endif
  const int min = 1, max = N / (THREAD_COUNT * 2);
  const int ignore = 0;
  const int n_opt = 4;
  const int n_iter = 20;

  /*PATSMA INSTANCE*/
  Autotuning *at = new Autotuning(min, max, ignore, dim, n_opt, n_iter);

  int *chunk = new int[dim];
  at->entireExecRuntime(matrix_calculation, chunk, A, N - 1);
  initialize(A, n + 1);

  printf("\n\n-----------------------Parallel Red Black Solver-----------------------\n\n\n");
  int iters;
  int convergence = false;
  double tmp;
  double diff;
  int i, j;
  for (iters = 1; iters < max_iterations; ++iters) {
    diff = 0;
    diff = matrix_calculation(A, N - 1, chunk);
    printf("Difference after %3d iterations: %f\n", iters, diff);
    if (diff / ((double)N * (double)N) < tolerance) {
      printf("\nConvergence achieved after %d iterations....Now exiting\n\n", iters);
      printf("PATSMA Points: %d %d\n", chunk[0], chunk[1]);
      at->print();
      return;
    }
  }
  printf("\n\nIteration LIMIT Reached...Exiting\n\n");
  printf("PATSMA Points: %d %d\n", chunk[0], chunk[1]);
  at->print();
}

int main(int argc, char *argv[]) {

  int i;
  long t_start, t_end;
  double time;
  A = new double *[N + 2];
  for (i = 0; i < N + 2; i++) {
    A[i] = new double[N + 2];
  }

  initialize(A, N);

  t_start = usecs();
  solve_parallel(A, N - 1);
  t_end = usecs();

  time = ((double)(t_end - t_start)) / 1000000;
  printf("Computation time for Parallel approach(secs): %f\n\n", time);
}