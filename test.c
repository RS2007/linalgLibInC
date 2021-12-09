#include <stdio.h>
#include "matrix.h"

int main(int argc, char *argv[])
{
	nml_mat *new = nml_mat_iden(3);
    nml_mat_print(new);
	return 0;
}
