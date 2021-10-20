#include <stdio.h>
#include "matrix.h"



int main(int argc, char *argv[])
{
    nml_mat* temp = nml_mat_iden(3);
    nml_mat* temp2 = nml_mat_iden(3);
    printf("%d ",nml_mat_equal(temp,temp2));
	return 0;
}
