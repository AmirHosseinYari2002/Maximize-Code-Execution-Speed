#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#pragma GCC optimize ("Ofast")

double fx(double x){
    double y = x;
    double ex = 1.0;

    for(int i=0; i<50; i++){
        ex += y;
        y  *= (x / (i+2));
    }
    return(ex);
}

double fx1(double x){
    register double y = 1.0;
    register double ex = 1.0 + x;
    register double k = 6.0;
    for(int i=0; i<25; i+=2){
        y *= x*x;
        ex += (y*x+(i+3)*y)/k;
        k *= (i+4)*(i+5);
    }
    return(ex);
}

double fx2(double x){
    register double y1 = x;
    register double ex = 1.0;
    register double y2, y3, y4;
    for(int i=0; i<50; i+=4){
        ex += y1;
        y2  = y1*(x / (i+2));
        ex += y2;
        y3  = y2*(x / (i+3));
        ex += y3;
        y4  = y3*(x / (i+4));
        ex += y4;
        y1  = y4*(x / (i+5));
    }
    return(ex);
}

double fx3(double* A,double* B, long int start,  long int count, long long int N, double K11, double K22){
    register double ex = 0.0;
    register long int iter;
    register double A_x, B_x;
    for(iter=start; iter <=count+start-1; iter++){
        // A
        A_x = iter * K11;
        __m256d A_y  = _mm256_set_pd(1, A_x, A_x*A_x/2, A_x*A_x*A_x/6);
        __m256d A_ex = _mm256_setzero_pd();
        // B
        B_x = iter * K22;
        __m256d B_y  = _mm256_set_pd(1,B_x,B_x*B_x/2,B_x*B_x*B_x/6);
        __m256d B_ex = _mm256_setzero_pd();
        for(register int i=0; i<50; i+=4){
            A_ex = _mm256_add_pd(A_y, A_ex);
            A_y  = _mm256_mul_pd(_mm256_set_pd((A_x*A_x*A_x*A_x)/((i+1)*(i+2)*(i+3)*(i+4)), (A_x*A_x*A_x*A_x)/((i+2)*(i+3)*(i+4)*(i+5)), (A_x*A_x*A_x*A_x)/((i+3)*(i+4)*(i+5)*(i+6)), (A_x*A_x*A_x*A_x)/((i+4)*(i+5)*(i+6)*(i+7))), A_y);
            B_ex = _mm256_add_pd(B_y, B_ex);
            B_y  = _mm256_mul_pd(_mm256_set_pd((B_x*B_x*B_x*B_x)/((i+1)*(i+2)*(i+3)*(i+4)), (B_x*B_x*B_x*B_x)/((i+2)*(i+3)*(i+4)*(i+5)), (B_x*B_x*B_x*B_x)/((i+3)*(i+4)*(i+5)*(i+6)), (B_x*B_x*B_x*B_x)/((i+4)*(i+5)*(i+6)*(i+7))), B_y);
        }
        // A * B
        ex += (A_ex[0]+A_ex[1]+A_ex[2]+A_ex[3]) * (B_ex[0]+B_ex[1]+B_ex[2]+B_ex[3]);
    }
    return (ex);
}

#define N 0x10000000

// main for fx3
int main(){

    clock_t t0, t1;
    int NT = 32;
    double *A  = (double*) malloc(N * sizeof(double));
    double *B  = (double*) malloc(N * sizeof(double));

    if(A == NULL || B == NULL){
        printf("Memory Allocation Error\n\n");
        return(-1);
    }

    register double sum = 0.0;
    long unsigned int i;

    const register double K1 = 3.0;
    const register double K2 = 2.0;
    register double K11 = K1 / N;
    register double K22 = K2 / N;

    t0 = clock();

    double E[NT][8];
    #pragma omp parallel num_threads(NT)
    {
        int id = omp_get_thread_num();
        E[id][8] += fx3(A, B, N/NT*id, N/NT, N, K11, K22);
    }

    for(int i=0; i<NT; i++){
        sum += E[i][8];
    }

    t1 = clock();

    printf("sum: %0.7f in %0.2f secs\n\n", sum,(float)(t1-t0)/CLOCKS_PER_SEC);
}

// // main for fx1 and fx2
//int main(){
//
//    clock_t t0, t1;
//    int NT = 8;
//    double *A  = (double*) malloc(N * sizeof(double));
//    double *B  = (double*) malloc(N * sizeof(double));
//
//    if(A == NULL || B == NULL){
//        printf("Memory Allocation Error\n\n");
//        return(-1);
//    }
//
//    register double sum = 0.0;
//    long unsigned int i;
//
//    t0 = clock();
//
//    const register double K1 = 3.0;
//    const register double K2 = 2.0;
//    const register double K11 = K1 / N;
//    const register double K22 = K2 / N;
//
//    #pragma omp parallel num_threads(NT)
//    {
//        long int i;
//        #pragma omp for reduction (+:sum)
//        for(i=0; i<N; i+=4){
//            __m256d vec_c = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
//            __m256d vec_a = _mm256_set_pd(fx1(i * K11), fx1((i+1) * K11), fx1((i+2) * K11), fx1((i+3) * K11));
//            __m256d vec_b = _mm256_set_pd(fx1(i * K22), fx1((i+1) * K22), fx1((i+2) * K22), fx1((i+3) * K22));
//            vec_c = _mm256_add_pd(_mm256_mul_pd(vec_a, vec_b), vec_c);
//            sum += vec_c[0] + vec_c[1] + vec_c[2] + vec_c[3] + vec_c[4];
//        }
//    }
//    t1 = clock();
//    printf("sum: %0.7f in %0.2f secs\n\n", sum,(float)(t1-t0)/CLOCKS_PER_SEC);
//}
