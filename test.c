#include <stdio.h>
#include "matrix.h"

int main(int argc, char *argv[])
{
    nml_mat *new = nml_mat_new(4, 4);
    new->data[0][0] = 1;
    new->data[0][1] = 2;
    new->data[0][2] = 2;
    new->data[0][3] = 2;
    new->data[1][0] = 2;
    new->data[1][1] = 5;
    new->data[1][2] = 2;
    new->data[1][3] = 5;
    new->data[2][0] = 2;
    new->data[2][1] = 2;
    new->data[2][2] = 9;
    new->data[2][3] = 1;
    new->data[3][0] = 2;
    new->data[3][1] = 5;
    new->data[3][2] = 1;
    new->data[3][3] = 7;
    nml_mat_lup *a = nml_mat_LU(new);
    nml_mat_print(a->L);
    nml_mat_print(a->U);
    nml_mat_print(a->P);

    return 0;
}
