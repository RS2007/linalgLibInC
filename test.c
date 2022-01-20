#include <stdio.h>
#include "matrix.h"

int main(int argc, char *argv[])
{
    nml_mat *A = nml_mat_new(2, 2);
    nml_mat *B = nml_mat_new(2, 1);
    A->data[0][0] = 2;
    A->data[0][1] = 5;

    A->data[1][0] = 4;
    A->data[1][1] = 5;


    B->data[0][0] = 20;
    B->data[1][0] = 10;

    nml_mat *x = solve_linear_LU(A, B);
//    nml_mat *x = nml_mat_mul_naive(A,B);
    nml_mat_print(x);
    nml_mat_free(x);
    nml_mat_free(A);
    nml_mat_free(B);
    return 0;
}
