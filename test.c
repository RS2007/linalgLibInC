#include <stdio.h>
#include "matrix.h"

int main(int argc, char *argv[])
{
	nml_mat *new = nml_mat_iden(3);
	new = nml_mat_multiply_scalar(new,5);
	nml_mat *new1 = nml_mat_iden(3);
    nml_mat_print(nml_mat_add(new,new1));
	return 0;
}
