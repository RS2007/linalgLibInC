#include <stdio.h>
#include "matrix.h"

int main(int argc, char *argv[])
{
    nml_mat *new = nml_mat_new(4, 5);
    new->data[0][0] = 3;
    new->data[0][1] = -2;
    new->data[0][2] = 5;
    new->data[0][3] = 0;
    new->data[0][4] = 2;
    new->data[1][0] = 4;
    new->data[1][1] = 5;
    new->data[1][2] = 8;
    new->data[1][3] = 1;
    new->data[1][4] = 4;
    new->data[2][0] = 1;
    new->data[2][1] = 1;
    new->data[2][2] = 2;
    new->data[2][3] = 1;
    new->data[2][4] = 5;
    new->data[3][0] = 2;
    new->data[3][1] = 7;
    new->data[3][2] = 6;
    new->data[3][3] = 5;
    new->data[3][4] = 7;
    nml_mat *ref = nml_mat_ref(new);
    nml_mat_print(ref);
    return 0;
}
