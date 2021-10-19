#include <stdio.h>
#include "matrix.h"



int main(int argc, char *argv[])
{
    FILE* fp;
    fp = fopen("./matrix01.txt","r");
    nml_mat* temp = nml_mat_fromfile(fp);
    int i,j;
    for(i=0;i<temp->num_rows;++i){
        for(j=0;j<temp->num_cols;++j){
            printf("%f ",temp->data[i][j]);
        }
    }
    fclose(fp);
	return 0;
}
