#include "task_28_11.h"

// распечатся матриц размера n
void print_matrix(double *A, int n) {
    printf("\n");
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%1.14f ", A[i * n + j]);
        }
        printf("\n");
    }
}

void generate_x_1_2(double a1, double a2, double *x1, double *x2, int i, double precision, int n){
    double norm;

    *x1 = a1 - sqrt(a1 * a1 + a2 * a2);
    *x2 = a2;

    norm = sqrt((*x1)*(*x1)+(*x2)*(*x2));

    // if(debug){
    //     printf("\na1 = %lf, a2 = %lf, x1 = %lf, x2 = %lf, norma = %lf.\n",a1,a2,x1,x2,norm);
    // }

    if (fabs(norm) <= precision){
        *x1 = 0;
        *x2 = 0;
    }
    else{
        *x1 /= norm;
        *x2 /= norm;                
    }
}

// перемножение матрицы U на матрицу А, в результате получится матрица, которая отличается только 2 строки от матрицы А
// поэтому выделяется и изменяется 2*n элементов. Важно, что в зависимости от шага будут меняться разные строки.
void multiply_U_A1(double *tmp, double *A, double u0, double u1, double u2, double u3, int i, int n, double precision){
	
	for (int j = i; j < n; j++) {
        tmp[n * n + j] = u0 * A[i * n + j] + u1 * A[i * n + n + j];
        tmp[n * n + n + j] = u2 * A[i * n + j] + u3 * A[i * n + n + j];
    }

    for (int j = i; j < n; j++) {
        A[i * n + j] = tmp[n * n + j];
        A[i * n + n + j] = tmp[n * n + n + j];

        if (fabs(A[i * n + n + j]) < precision)
            A[i * n + n + j] = 0;
    }

}

// произведение матриы Q на U.
void multiply_Q_U(double *tmp, double u0, double u1, double u2, double u3, int step, int n) {
    int i;
    for (i = 0; i < n; i++) {
        tmp[n * n + i] = tmp[i * n + step] * u0 + tmp[i * n + step + 1] * u1;
        tmp[n * n + n + i] = tmp[i * n + step] * u2 + tmp[i * n + step + 1] * u3;
    }
    for (i = 0; i < n; i++) {
        tmp[i * n + step] = tmp[n * n + i];
        tmp[i * n + step + 1] = tmp[n * n + n + i];
    }
}

// сортировка пузырьком значений матрицы E, где хранятся собственные значаения для матрицы А
void sort(double *E, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (E[i] < E[j]) {
                double temp = E[j];
                E[j] = E[i];
                E[i] = temp;
            }
        }
    }
}

// Умножение матрицы A ha Q, 
void multiply_A_Q(double *A, double *tmp, int n, double epsilon) {
    int i, j, k;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            tmp[n * n + 2 * n + i * n + j] = 0.0;
        }
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                tmp[n * n + 2 * n + i * n + j] += A[i * n + k] * tmp[k * n + j];
            }
        }
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A[i * n + j] = tmp[n * n + 2 * n + i * n + j];
        }
    }
}

// проверка на выход из цикла while - сравнение собственных значаений на предыдущем и нынешнем шаге с заданным параметром 
// epsilon. Если хотя бы один элемент больше на epsilon, то возвращает 1, иначе 0
int check_diff(double *A, double *E, int n, double epsilon) {
    int i;
    for (i = 0; i < n; i++) {
        if (fabs(A[i * n + i] - E[i]) > epsilon) {
            return 1;
        }
    }
    return 0;
}

// проверка на выход из цикла while - проверка матрицы Q
int check_Q(double *tmp, int n, double epsilon) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                if (fabs(tmp[i * n + i]) > epsilon)
                    return 1;
            }
        }
    }
    return 0;
}

// создается или обновляется матрица Q, которая является единичной - I
void mod_Q(double *tmp, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                tmp[i * n + j] = 0.0;
            } else {
                tmp[i * n + j] = 1.0;
            }
        }
    }
}

// обновляются значения вектора Е
void rewrite_E(double *A, double *E, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        E[i] = A[i * n + i];
    }
}

int evc_28_11(int n, int max_iterations, double epsilon, double *A, double *E, double *tmp, double precision) {

    int flag_d=0;    // флаг для выхода при достижении дигональной матрицы
    int flag_q=0;    // разность между собственными значениями матрицми меньше заданного epsilon, и можно завершить алгоритм по условию
    int number_iter = 0; // количество итераций
    int status = 0;    //служит как флаш для пропуска шага цикла
    int i, j;

    double x1, x2;
    mod_Q(tmp, n);

    double u0, u1, u2, u3;

    while (max_iterations == 0 || number_iter < max_iterations) {

        flag_d = 0;
        flag_q = 0;

        for (int i = 0; i < n - 1; i++) {
            // если элемент 0, то преобразоваение не нужно делать
            if (fabs(A[i * n + n + i]) <= precision){
                A[i * n + n + i] = 0.0;
            }
            else{
                ////////////////////////////////////////////////
                generate_x_1_2(A[i * n + i], A[i * n + n + i] , &x1, &x2, i, precision, n);
                ////////////////////////////////////////////////
    // определяется матрца U
                u0 = 1 - 2 * x1 * x1;
                u3 = 1 - 2 * x2 * x2;
                u1 = -2 * x1 * x2;
                u2 = u1;

                ///////////////////////////
                multiply_U_A1( tmp, A,  u0,  u1,  u2,  u3,  i,  n,  precision);
                //////////////////////////
                multiply_Q_U(tmp, u0, u1, u2, u3, i, n);
            }
        }

        if (debug == 1) {
            printf("\n");
            printf("Evc ");
            printf("iteration number = %d", number_iter);

            printf("\n");
            printf("R matrix");
            print_matrix(A, n);

            printf("\n");
            printf("Q matrix");
            print_matrix(tmp, n);
        }

        multiply_A_Q(A, tmp, n, epsilon);

        flag_q = check_Q(tmp, n, epsilon);
        flag_d = check_diff(A, E, n, epsilon);

        mod_Q(tmp, n);
        rewrite_E(A, E, n);
// выполняется одно из условий выхода из цикла while
        if (flag_d * flag_q == 0) {
            break;
        }
        number_iter++;
    }

    sort(E, n);
// при условии, что в матрице Q есть отличия в элементах вне диагонали и нужной точности 
    //в i и i+1 шагах не достигнута
    if (flag_d * flag_q == 1)
        return 1;

    return 0;
}


int evc_memsize_28_11(int n) {
    return n * n + n * n + 2 * n;
}

