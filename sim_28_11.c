#include "task_28_11.h"

//// Строим вектор х, считая норму вектора х в двух переменных уменьшаем 
//вычисления так как разница в нормах a_i и a_i-e_1||a_i|| в одном элементе при вычислении
void generate_x (double *A, double *tmp,int i,int n){
    double norma_part=0;
    double norm = 0;
    double prec =1e-14;
    int debug=0;
    int j,k; 

    for(j=0;j<n-i-1;j++){
        tmp[j]=A[i*n+j+i+1];
        if(j!=0){
            norma_part+=tmp[j]*tmp[j];
        }
    }
    tmp[0]-=sqrt(norma_part+tmp[0]*tmp[0]);
    norm = sqrt(norma_part+tmp[0]*tmp[0]);
    
    for( k=0;k<n-i-1;k++){
        if( norm < prec){
            tmp[k]=0;
        }
        else{
            tmp[k]/=norm;
        }
    }
    if(debug){
        printf("Generating norm vector for step %d \n",i);
        for( k=0;k<n-i-1;k++){
            printf("%1.9lf ",tmp[k] );
        }
        printf("\n");
        printf("Norm equal to = %1.9lf \n",norm);
    }
}
// генерируем матрицу U, как U=I-2x*x_t, зная что она отличается от I только в квадратной матрице
void generate_U(double *tmp,int i, int n){
    int j,k;
	for(j=0;j<n;j++){
        for( k=0;k<n;k++){
            if(j==k){
                if(j<i+1){
                    tmp[n+j*n+k]=1;
                }
                else{
                    tmp[n+j*n+k]=1-2*tmp[j-(i+1)]*tmp[j-(i+1)];
                }
            }
            else{
                if(j<i+1 || k<i+1){
                    tmp[n+j*n+k]=0;
                }
                else{
                    tmp[n+j*n+k]=-2*tmp[j-(i+1)]*tmp[k-(i+1)];
                }
            }
        }
    }
}
// копирование матрицу А в tmp, матрицу А зануляем
void copy_matrix_A_tmp(double *A,double *tmp,int n){
    int index=n+n*n;
	int i,j;
    for( i=0;i<n;i++){
        for(j=0;j<n;j++){
            tmp[index+i*n+j]=A[i*n+j];
            A[i*n+j]=0;
        }
    }
    
}
// умножение матрицы U на А, результат записывается в А
void multiply_U_A(double *tmp, double*A,int n){
	int i,j,k;
   for( i=0;i<n;i++){
       for( j=0;j<n;j++){
           for( k=0;k<n;k++){
               A[n*i+j]+=tmp[n+i*n+k]*tmp[n+n*n+k*n+j];
           }
       }
   }
   
}
// умножение матрицы A_1=U*A на  U и получается A_k+1, результат записывается в А
void multiply_UA_U(double *A, double*tmp,int n){
	int i,j,k;
    for( i=0;i<n;i++){
        for( j=0;j<n;j++){
            tmp[n+n*n+n*i+j]=0;
        }
    }
    
   for( i=0;i<n;i++){
       for( j=0;j<n;j++){
           for( k=0;k<n;k++){
               tmp[n+n*n+n*i+j]+=A[i*n+k]*tmp[n+k*n+j];
           }
       }
   }
   
   for( i=0;i<n;i++){
       for( j=0;j<n;j++){
           A[i*n+j]=tmp[n+n*n+n*i+j];
       }
   }
}

int sim_28_11(int n, double* A, double* tmp, double precision){  
	int i;
    for( i=0;i<n-2;i++){
        generate_x(A,tmp,i,n);
        generate_U(tmp,i,n);
        copy_matrix_A_tmp(A,tmp,n);
        multiply_U_A(tmp,A,n);
        multiply_UA_U(A,tmp,n);
    }
    
    return 0;
}

int sim_memsize_28_11(int n){
    return n+n*n+n*n;
}
