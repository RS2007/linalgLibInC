#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "matrix.h"
#include <string.h>

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
	unsigned int i;

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
	unsigned int i;
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

double nml_rand_interval(int min, int max)
{
	double d = (double)((rand() % (max - min + 1)) + min);
	return d;
}

nml_mat *nml_mat_rnd(unsigned int num_rows, unsigned int num_cols, int min, int max)
{
	nml_mat *res = nml_mat_new(num_rows, num_cols);
	unsigned int i, j;

	for (i = 0; i < num_rows; ++i)
	{
		for (j = 0; j < num_cols; ++j)
		{
			res->data[i][j] = nml_rand_interval(min, max);
		}
	}
	return res;
}

nml_mat *nml_mat_sqr(unsigned int rowCol)
{
	return nml_mat_new(rowCol, rowCol);
}

nml_mat *nml_mat_sqr_rnd(unsigned int rowCol, int min, int max)
{
	return nml_mat_rnd(rowCol, rowCol, min, max);
}

nml_mat *nml_mat_iden(unsigned int rowCol)
{
	nml_mat *iden = nml_mat_sqr(rowCol);
	int i;
	for (i = 0; i < rowCol; ++i)
	{
		iden->data[i][i] = 1.0;
	}
	return iden;
}

nml_mat *nml_mat_fromfile(FILE *f)
{
	int i, j;
	unsigned int num_rows = 0;
	unsigned int num_cols = 0;
	fscanf(f, "%u", &num_rows);
	fscanf(f, "%u", &num_cols);
	nml_mat *res = nml_mat_new(num_rows, num_cols);
	for (i = 0; i < res->num_rows; ++i)
	{
		for (j = 0; j < res->num_cols; ++j)
		{
			fscanf(f, "%lf\t", &res->data[i][j]);
		}
	}
	return res;
}
int nml_mat_equaldim(nml_mat *m1, nml_mat *m2)
{
	if (m1->num_cols == m2->num_cols && m1->num_rows == m2->num_rows)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

int nml_mat_equal(nml_mat *m1, nml_mat *m2)
{
	if (nml_mat_equaldim(m1, m2) == 0)
	{
		printf("Dimensions of matrix unequal");
		return 0;
	}
	else
	{
		int i, j;
		int flag = 1;
		for (i = 0; i < m1->num_rows; ++i)
		{
			for (j = 0; j < m1->num_rows; ++j)
			{
				if (m1->data[i][j] == m2->data[i][j])
				{
					continue;
				}
				else
				{
					flag = 0;
					break;
				}
			}
		}
		return flag;
	}
}

nml_mat *nml_mat_get_column(nml_mat *matrix, unsigned int col)
{
	int i;
	nml_mat *column = nml_mat_new(matrix->num_rows, 1);
	for (i = 0; i < matrix->num_rows; ++i)
	{
		column->data[i][0] = matrix->data[i][col - 1];
	}
	return column;
}

nml_mat *nml_mat_get_row(nml_mat *matrix, unsigned int row)
{
	nml_mat *rowMat = nml_mat_new(1, matrix->num_cols);
	//row is contigous set of memory
	//syntax of memcpy is void* mecmcpy(void* to,const void* from,size_t )
	memcpy(rowMat->data[0], matrix->data[row - 1], matrix->num_cols * sizeof(double));
	return rowMat;
}

nml_mat *nml_set_all_elements(nml_mat *matrix, double num)
{
	int i, j;
	nml_mat *copy = matrix;
	for (i = 0; i < matrix->num_rows; ++i)
	{
		for (j = 0; j < matrix->num_cols; ++j)
		{
			copy->data[i][j] = num;
		}
	}
	return copy;
}

nml_mat *nml_set_diagonal_elements(nml_mat *matrix, double num)
{
	int i;
	nml_mat *copy = matrix;
	for (i = 0; i < matrix->num_rows; ++i)
	{
		copy->data[i][i] = num;
	}
	return copy;
}

nml_mat *nml_row_multipy_scalar(nml_mat *matrix, unsigned int row, int scalar)
{
	int i;
	for (i = 0; i < matrix->num_cols; ++i)
	{
		matrix->data[row][i] *= scalar;
	}
	return matrix;
}
nml_mat *nml_col_multiply_scalar(nml_mat *matrix, unsigned int col, int scalar)
{
	int i;
	for (i = 0; i < matrix->num_cols; ++i)
	{
		matrix->data[i][col] *= scalar;
	}
	return matrix;
}

nml_mat *nml_rows_add(nml_mat *matrix, unsigned int addendum, unsigned int original, int multiplier)
{
	int i;
	nml_mat *copy = matrix;
	for (i = 0; i < matrix->num_cols; ++i)
	{
		copy->data[original - 1][i] = copy->data[original - 1][i] + multiplier * copy->data[addendum - 1][i];
	}
	return copy;
}

void nml_mat_print(nml_mat *matrix)
{
	int i, j;
	for (i = 0; i < matrix->num_rows; ++i)
	{
		for (j = 0; j < matrix->num_cols; ++j)
		{
			printf(" %lf ", matrix->data[i][j]);
		}
		printf("\n");
	}
}

nml_mat *nml_mat_multiply_scalar(nml_mat *matrix, int scalar)
{
	nml_mat *copy = matrix;
	int i, j;
	for (i = 0; i < matrix->num_rows; ++i)
	{
		for (j = 0; j < matrix->num_cols; ++j)
		{
			copy->data[i][j] = copy->data[i][j] * scalar;
		}
	}
	return copy;
}

nml_mat *nml_mat_remove_column(nml_mat *matrix, unsigned int column)
{
	if (column > matrix->num_cols)
	{
		perror("Invalid columns");
		return NULL;
	}
	nml_mat *r = nml_mat_new(matrix->num_rows, matrix->num_cols - 1);
	int i, j, k;
	for (i = 0; i < matrix->num_rows; ++i)
	{
		for (j = 0, k = 0; j < matrix->num_cols; ++j)
		{
			if (column != j)
			{
                // k resets back to 0 after i++ [one j is missing and that's going to lead to an incomplete array]
                // but k has to be contiguous
				r->data[i][k++] = matrix->data[i][j];
			}
		}
	}
	return r;
}

nml_mat *nml_mat_remove_row(nml_mat* matrix,unsigned int row){
    if(row > matrix->num_rows){
        perror("Invalid rows");
        return NULL;
    }
    nml_mat* r = nml_mat_new(matrix->num_rows-1,matrix->num_cols);
    int i,j,k;
    for(i = 0,k = 0;i<matrix->num_rows;++i){
       for(j = 0;j<matrix->num_cols;++j){
           if(row != i){
               r->data[k][j] = matrix->data[i][j];
           }
       }
       k++;
    }
    return r;
}

nml_mat* nml_mat_add(nml_mat* matrix1,nml_mat* matrix2){
	if(matrix1->num_rows != matrix2->num_rows || matrix1->num_cols != matrix2->num_cols){
	perror("Cannot add matrices of different dimensions");
	return NULL;	
	}
	nml_mat* sum = nml_mat_new(matrix1->num_rows,matrix1->num_cols);
	int i,j;
	for(i = 0;i<matrix1->num_rows;++i){
		for(j = 0;j<matrix1->num_cols;++j){
			sum->data[i][j] = matrix1->data[i][j] + matrix2->data[i][j];	
		}
	}	
	return sum;
}
	


