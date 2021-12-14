#include <stdio.h>
#include "matrix.h"

int main(int argc, char *argv[])
{
    nml_mat *new = nml_mat_new(3, 3);
    new = nml_set_all_elements(new, 5);
    int k = nml_mat_get_col_pivot(new, 1);
    printf("%d", k);
    return 0;
}
