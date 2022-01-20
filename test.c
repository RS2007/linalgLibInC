#include <stdio.h>
#include "matrix.h"

int main(int argc, char *argv[])
{
    nml_mat *A = nml_mat_new(2, 2);
    nml_mat *B = nml_mat_new(2, 1);
    A->data[0][0] = 1;
    A->data[0][1] = 1;

    A->data[1][0] = 2;
    A->data[1][1] = 3;

    B->data[0][0] = 13;
    B->data[1][0] = -1;

    nml_mat *x = solve_linear_LU(A, B);
    nml_mat_print(x);
    return 0;
}