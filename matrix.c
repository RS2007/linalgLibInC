#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "matrix.h"


nml_mat *nml_mat_new(unsigned int num_rows, unsigned int num_cols)
{
	//Check if rows are 0
	if (num_rows == 0)
	{
		printf("Invalid rows");
		return NULL;
	}

	//Check if columns are 0
	if (num_cols == 0)
	{
		printf("Invalid columns");
		return NULL;
	}

	//giving nml_mat its memory
	nml_mat *m = (nml_mat *)malloc(sizeof(nml_mat));
	//NP_CHECK(m);
	m->num_rows = num_rows;
	m->num_cols = num_cols;
	m->isSquare = (num_rows == num_cols) ? 1 : 0;

	//data is pointer to pointer to double(set of arrays which are a set of doubles)
	//setting the outer layer(row layer)

	m->data = calloc(m->num_rows, sizeof(*m->data));
	//NP_CHECK(m->data);
	int i;

	//structure established complete structure by adding the columns/values

	for (i = 0; i < m->num_rows; ++i)
	{

		//double de-referencing for individual numbers iterating over rows

		m->data[i] = calloc(m->num_cols, sizeof(**m->data));
		//NP_CHECK(m->data[i]);
	}
	return m;
}
void nml_mat_free(nml_mat *matrix)
{
	int i;
	//each rows memory is freed
	for (i = 0; i < matrix->num_rows; ++i)
	{
		free(matrix->data[i]);
	}
	//freeing the outer structure
	free(matrix->data);
	//freeing the whole matrix
	free(matrix);
}


double nml_rand_interval(int min,int max){
    double d = (double)((rand() % (max-min+1))+min);
	return d;
}


nml_mat* nml_mat_rnd(unsigned int num_rows,unsigned int num_cols,int min,int max){
	nml_mat* res = nml_mat_new(num_rows,num_cols);
	int i,j;
	
	for(i=0;i<num_rows;++i){
		for(j=0;j<num_cols;++j){
			res->data[i][j] = nml_rand_interval(min,max);
		}	
	}	
	return res;
}

nml_mat* nml_mat_sqr(unsigned int rowCol){
    return nml_mat_new(rowCol,rowCol);
}

nml_mat* nml_mat_sqr_rnd(unsigned int rowCol,int min,int max){
    return nml_mat_rnd(rowCol,rowCol,min,max);
}

nml_mat* nml_mat_iden(unsigned int rowCol){
    nml_mat* iden = nml_mat_sqr(rowCol);
    int i;
    for(i = 0; i<rowCol;++i){
        iden->data[i][i] = 1.0;
    }
    return iden;
}

nml_mat* nml_mat_fromfile(FILE *f){
    int i,j;
    unsigned int num_rows=0;
    unsigned int num_cols=0;
    fscanf(f,"%u",&num_rows);
    fscanf(f,"%u",&num_cols);
    nml_mat *res = nml_mat_new(num_rows,num_cols);
    for(i = 0;i<res->num_rows;++i){
        for(j=0;j<res->num_cols;++j){
            fscanf(f,"%lf\t",&res->data[i][j]);
        }
    }
    return res;
}
int nml_mat_equaldim(nml_mat* m1,nml_mat* m2){
    if(m1->num_cols == m2->num_cols && m1->num_rows == m2->num_rows){
        return 1;
    }else{
        return 0;
    }
}

int nml_mat_equal(nml_mat* m1,nml_mat* m2){
    if(nml_mat_equaldim(m1,m2) == 0){
        printf("Dimensions of matrix unequal");
        return 0;
    }else{
       int i,j;
       int flag=1;
       for(i=0;i<m1->num_rows;++i){
           for(j=0;j<m1->num_rows;++j){
               if(m1->data[i][j] == m2->data[i][j]){
                   continue;
               }else{
                   flag = 0;
                   break;
               }
           }
       }
       return flag;
    }
}
