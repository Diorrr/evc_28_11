#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int debug, error;
int evc_memsize_28_11(int n);
int sim_memsize_28_11(int n);

int sim_28_11(int n, double* A, double* tmp, double precision);
int evc_28_11(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision);

