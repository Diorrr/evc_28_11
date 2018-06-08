
#include "task_28_11.h"

int main(int argc, char const *argv[]) {
    const char *inFileName = "28_11_in.txt";
    const char *outFileName = "28_11_out.txt";

    double precision = 1e-14;
    int max_iterations = 0;
    double epsilon = 1e-10;

    //error = 0;
    int p = 0, t = 0, h = 0;

    double *A;
    double *tmp;
    double *tmp1;
    double *E;

    int i = 1;
    int j;
    int n;
    int status;

    if (argc > i) {
        if (argv[i][0] != '-') {
            inFileName = argv[i++];
        }
    }
    if (argc > i) {
        if (argv[i][0] != '-') {
            outFileName = argv[i++];
        }
    }

    int in_file_read = 0;
    int out_file_read = 0;

    for (; i < argc; i++) {
        char const *arg = argv[i];

        if (arg[0] == '-' && arg[2] == '\0') {

            if (arg[1] == 'd') {
                debug = 1;
                continue;
            } else if (arg[1] == 'e') {
                error = 1;
                continue;
            } else if (arg[1] == 'p') {
                p = 1;
                continue;
            } else if (arg[1] == 't') {
                t = 1;
                continue;
            } else if (arg[1] == 'h' || arg[1] == '?') {
                h = 1;
                continue;
            } else {
                return 1;
            }
        }

        if (sscanf(arg, "-prec=%lf%*c", &precision) == 1) {
            continue;
        } else if (sscanf(arg, "-eps=%lf%*c", &epsilon) == 1) {
            continue;
        } else if (sscanf(arg, "-max_iter=%d%*c", &max_iterations) == 1) {
            continue;
        } else if (argv[i][0] != '-' && in_file_read == 0) {
            inFileName = argv[i];
            in_file_read = 1;
        } else if (argv[i][0] != '-' && out_file_read == 0) {
            outFileName = argv[i];
            out_file_read = 1;
        } else {
            return 1;
        }
    }

    if (h) {
        printf("Usage: evc [input_file_name] [output_file_name] [options]\n"
               "Where options include:\n"
               "  -d    print debug messages [default OFF]\n"
               "  -e    print errors [default OFF]\n"
               "  -p    print matrix [default OFF]\n"
               "  -t    print execution time [default OFF]\n"
               "  -prec=<num>       precision [default - 1e-14]\n"
               "  -eps=<num>        'epsilon' [default - 1e-10]\n"
               "  -max_iterations=<num>   limited number of iterations\n"
               "    [default - 0, i.e. not limit]\n"
               "    -h, -?     print this and exit\n");
        return 0;
    }


    FILE *file_in = fopen(inFileName, "r");
    FILE *file_out = fopen(outFileName, "w");

    if (file_in == NULL) {
        printf("\n can not open file");
        return 0;
    }

    fscanf(file_in,"%d", &n);

    A = (double *) malloc((n * n) * sizeof(double));
    tmp = (double *) malloc(sim_memsize_28_11(n) * sizeof(double));
    tmp1 = (double *) malloc(evc_memsize_28_11(n) * sizeof(double));
    E = (double *) malloc(n * sizeof(double));

    for ( i = 0;i < n;i++) {
        for (j = 0;j < n;j++) {
            fscanf(file_in,"%lf", &A[i*n+ j]);
        }
    }

    if (p) {
        for (i = 0;i < n;i++) {
            for (j = 0; j < n; j++) {
                printf("% 1.9lf ", A[i * n+ j]);
            }
            printf("\n");
        }
    }

    clock_t time_c = clock();

    status = sim_28_11(n, A, tmp, precision);
//// Result of sim
    if (debug) {
        printf("\n");
        for (i = 0;i < n;i++) {
            for (j = 0;j < n;j++) {
                printf("% 1.9lf ", A[i*n+ j]);
            }
            printf("\n");
        }
    }
////Result of evc
    status = evc_28_11(n, max_iterations, epsilon, A, E, tmp1, precision);
    if (status != 0) {
        fprintf(file_out,
                "0");
        printf("\n Not Success");
        return status;
    }
    else if (status == 0) {
        fprintf(file_out,"%d\n", n);
        for (i = 0; i < n;i++) {
            fprintf(file_out,"%1.9lf\n", E[i]);
        }
        printf("\n Success");

    }
    if (t) {
        printf("time: %lf seconds\n", (double) (clock()- time_c) / CLOCKS_PER_SEC);
    }

//    fprintf(file_out, "%d\n", n);

//    for (i = 0; i < n; i++) {
//        fprintf(file_out, "%1.9lf\n", E[i]);
//    }

    free(A);
    free(tmp);
    free(E);

    fclose(file_in);
    fclose(file_out);

    return 0;

}
