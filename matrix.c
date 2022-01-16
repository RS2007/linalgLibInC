#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include "matrix.h"
#include <string.h>

#define NML_MIN_COEFF 0.00001

nml_mat *nml_mat_new(unsigned int num_rows, unsigned int num_cols)
{
	// Check if rows are 0
	if (num_rows == 0)
	{
		printf("Invalid rows");
		return NULL;
	}

	// Check if columns are 0
	if (num_cols == 0)
	{
		printf("InvaliGd columns");
		return NULL;
	}

	// giving nml_mat its memory
	nml_mat *m = (nml_mat *)malloc(sizeof(nml_mat));
	// NP_CHECK(m);
	m->num_rows = num_rows;
	m->num_cols = num_cols;
	m->isSquare = (num_rows == num_cols) ? 1 : 0;

	// data is pointer to pointer to double(set of arrays which are a set of doubles)
	// setting the outer layer(row layer)

	m->data = calloc(m->num_rows, sizeof(*m->data));
	// NP_CHECK(m->data);
	unsigned int i;

	// structure established complete structure by adding the columns/values

	for (i = 0; i < m->num_rows; ++i)
	{

		// double de-referencing for individual numbers iterating over rows

		m->data[i] = calloc(m->num_cols, sizeof(**m->data));
		// NP_CHECK(m->data[i]);
	}
	return m;
}
void nml_mat_free(nml_mat *matrix)
{
	unsigned int i;
	// each rows memory is freed
	for (i = 0; i < matrix->num_rows; ++i)
	{
		free(matrix->data[i]);
	}
	// freeing the outer structure
	free(matrix->data);
	// freeing the whole matrix
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
	// row is contigous set of memory
	// syntax of memcpy is void* mecmcpy(void* to,const void* from,size_t )
	memcpy(rowMat->data[0], matrix->data[row - 1], matrix->num_cols * sizeof(double));
	return rowMat;
}

nml_mat *nml_set_all_elements(nml_mat *matrix, double num)
{
	int i, j;
	nml_mat *copy = nml_mat_cp(matrix);
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
	nml_mat *copy = nml_mat_cp(matrix);
	for (i = 0; i < matrix->num_rows; ++i)
	{
		copy->data[i][i] = num;
	}
	return copy;
}

nml_mat *nml_row_multipy_scalar(nml_mat *matrix, unsigned int row, double scalar)
{
	int i;
	for (i = 0; i < matrix->num_cols; ++i)
	{
		matrix->data[row][i] *= scalar;
	}
	return matrix;
}
nml_mat *nml_col_multiply_scalar(nml_mat *matrix, unsigned int col, double scalar)
{
	int i;
	for (i = 0; i < matrix->num_cols; ++i)
	{
		matrix->data[i][col] *= scalar;
	}
	return matrix;
}

nml_mat *nml_rows_add(nml_mat *matrix, unsigned int addendum, unsigned int original, double multiplier)
{
	int i;
	nml_mat *copy = nml_mat_cp(matrix);
	for (i = 0; i < matrix->num_cols; ++i)
	{
		copy->data[original][i] = copy->data[original][i] + multiplier * copy->data[addendum][i];
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
	printf("\n\n");
}
nml_mat *nml_mat_multiply_scalar(nml_mat *matrix, int scalar)
{
	nml_mat *copy = nml_mat_cp(matrix);
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

nml_mat *nml_mat_remove_row(nml_mat *matrix, unsigned int row)
{
	if (row > matrix->num_rows)
	{
		perror("Invalid rows");
		return NULL;
	}
	nml_mat *r = nml_mat_new(matrix->num_rows - 1, matrix->num_cols);
	int i, j, k;
	for (i = 0, k = 0; i < matrix->num_rows; ++i)
	{
		for (j = 0; j < matrix->num_cols; ++j)
		{
			if (row != i)
			{
				r->data[k][j] = matrix->data[i][j];
			}
		}
		k++;
	}
	return r;
}

nml_mat *nml_mat_swap_row(nml_mat *matrix, unsigned int row1, unsigned int row2)
{
	if (row1 > matrix->num_rows || row2 > matrix->num_rows)
	{
		perror("Invalid row");
	}
	nml_mat *copy = nml_mat_cp(matrix);
	// swap the rows(contiguous memory locations)
	double *temp = copy->data[row1];
	copy->data[row1] = copy->data[row2];
	copy->data[row2] = temp;
	return copy;
}

nml_mat *nml_mat_swap_column(nml_mat *matrix, unsigned int col1, unsigned int col2)
{
	if (col1 > matrix->num_cols || col2 > matrix->num_cols)
	{
		perror("Invalid column");
	}

	nml_mat *copy = nml_mat_cp(matrix);
	// swap the numbers between the rows
	int i;
	for (i = 0; i < copy->num_rows; ++i)
	{
		double temp = copy->data[i][col1];
		copy->data[i][col1] = copy->data[i][col2];
		copy->data[i][col2] = temp;
	}

	return copy;
}

nml_mat *nml_mat_cp(nml_mat *matrix)
{
	int i, j;
	nml_mat *copy = nml_mat_new(matrix->num_rows, matrix->num_cols);
	for (i = 0; i < matrix->num_rows; ++i)
	{
		for (j = 0; j < matrix->num_cols; ++j)
		{
			copy->data[i][j] = matrix->data[i][j];
		}
	}
	return copy;
}

nml_mat *nml_mat_add(nml_mat *matrix1, nml_mat *matrix2)
{
	if (matrix1->num_rows != matrix2->num_rows || matrix1->num_cols != matrix2->num_cols)
	{
		perror("Cannot add matrices of different dimensions");
		return NULL;
	}
	nml_mat *sum = nml_mat_new(matrix1->num_rows, matrix1->num_cols);
	int i, j;
	for (i = 0; i < matrix1->num_rows; ++i)
	{
		for (j = 0; j < matrix1->num_cols; ++j)
		{
			sum->data[i][j] = matrix1->data[i][j] + matrix2->data[i][j];
		}
	}
	return sum;
}

nml_mat *nml_mat_subtract(nml_mat *matrix1, nml_mat *matrix2)
{
	if (matrix1->num_rows != matrix2->num_rows || matrix1->num_cols != matrix2->num_cols)
	{
		perror("Cannot subtract matrices of different dimensions");
		return NULL;
	}
	nml_mat *diff = nml_mat_new(matrix1->num_rows, matrix1->num_cols);
	int i, j;
	for (i = 0; i < matrix1->num_rows; ++i)
	{
		for (j = 0; j < matrix1->num_cols; ++j)
		{
			diff->data[i][j] = matrix1->data[i][j] - matrix1->data[i][j];
		}
	}
	return diff;
}

nml_mat *nml_mat_mul_naive(nml_mat *matrix1, nml_mat *matrix2)
{
	if (matrix1->num_cols != matrix2->num_rows)
	{
		perror("Cannot multiply matrices, incompatible dimensions");
		return NULL;
	}
	nml_mat *mul = nml_mat_new(matrix1->num_rows, matrix2->num_cols);
	int i, j, k;
	for (i = 0; i < matrix1->num_cols; ++i)
	{
		for (j = 0; j < matrix2->num_rows; ++j)
		{
			mul->data[i][j] = 0;
			for (k = 0; k < matrix1->num_cols; ++k)
			{
				mul->data[i][j] += matrix1->data[i][k] * matrix2->data[k][j];
			}
		}
	}
	return mul;
}

nml_mat *nml_mat_mul_strassen(nml_mat *matrix1, nml_mat *matrix2)
{
	if (matrix1->num_cols != matrix2->num_rows)
	{
		perror("Cannot multiply matrices, incompatible dimensions");
		return NULL;
	}
	return nml_mat_iden(3);
	// Filled by strassen algorithm
}

int nml_mat_get_col_pivot(nml_mat *matrix, unsigned int col, unsigned int row)
{
	int i;
	for (i = row; i < matrix->num_rows; ++i)
	{
		if (fabs(matrix->data[i][col]) > NML_MIN_COEFF)
		{
			return i;
		}
	}
	return -1;
}

nml_mat *nml_mat_ref(nml_mat *matrix)
{
	/* trying to make sense of what they are doing
		unsigned int currentRow = 0
		for( unsigned int c = 0;c<m->cols;++c){
			unsigned int r = currentRow;
			if(r >= m->rows){
				break; when the guys stupid(error handling)
			}
			for(;r<m->rows;r++){
				if(m->elements[r][c] != 0.0f){
					break
				}
			}
			if(r == m->rows){
				continue;
			}
			swapRows(m,currentRow+1,r+1); checking the first non zero and then swapping if necessary for right zeroes structure

			float factor = 1/m->elements[currentRow][c];
			for (unsigned int col = c;col<m->cols;++c){
				m->elements[currentRow][col] *= factor;
			}

			for(r = currentRow+1;r<m->rows;++r){
				addMultiple(m,r+1,currentRow+1,-1*m->elements[r][c]);
			}
			currentRow++;
		}
	*/

	// Row echelon form(staircase form)
	//  -- all zeroes at bottom
	//  -- first nonzero entry from left is a 1(leading 1)
	//  -- each leading 1 is to the right of the leading 1s in rows above

	// Gaussian algorithm - convert matrices to REF()
	//  -- If matrix all zeroes in REF
	//  -- For each column
	//  -- Find first row from top containing nonzero entry , move to top
	//  -- Multiplly new top row by 1/a to get leading 1
	//  -- subtract multiples of that row from rows below it to get zeros under leading 1
	nml_mat *copy = nml_mat_cp(matrix);
    int i,j;
    for(i = 0;i<copy->num_cols;++i){
        int pivot = nml_mat_get_col_pivot(copy,i,i);
        if(pivot<0)
            continue;
        if(pivot != i)
            copy = nml_mat_swap_row(copy,i,pivot);
        for(j = i+1;j<copy->num_rows;++j){
                copy = nml_rows_add(copy,i,j,-(copy->data[j][i])/copy->data[i][i]);
        }

    }
    return copy;
	// !REFACTOR change into a for loop
    
	//	num_rows,num_cols = arr.shape
	// if(num_rows != num_cols): return "error"
	// for k in range(num_rows-1):
	//		for i in range(k+1,num_rows):
	//			if arr[i,k] == 0: continue
	//			factor = arr[k,k]/arr[i,k]
	//			for j in range(k,num_rows):
	//				arr[i,j] = arr[k,j] - arr[i,j]*factor
}

nml_mat *nml_mat_rref(nml_mat *matrix)
{
	nml_mat *copy = nml_mat_cp(matrix);
	int i,j;
	for(i = 0;i<copy->num_cols;++i){
		int pivot = nml_mat_get_col_pivot(copy,i,i);
		if(pivot<0)
			continue;
		if(pivot != i)		
			copy = nml_mat_swap_row(copy,i,pivot);
		copy = nml_row_multipy_scalar(copy,i,1/copy->data[i][i]);
		nml_mat_print(copy);
		for(j = 0;j<copy->num_rows;++j){
			if(j != i){
				copy = nml_rows_add(copy,i,j,-(copy->data[j][i]));
			}	
			continue;
		}
	}
	return copy;
}

nml_mat_lup* nml_mat_lup_new(nml_mat* L,nml_mat* U,nml_mat* P,unsigned int num_permutations){
    nml_mat_lup *r = malloc(sizeof(*r));
    r->L = L;
    r->U = U;
    r->P = P;
    r->num_permutations= num_permutations;
    return r;
}

void nml_mat_lup_free(nml_mat_lup* LUP){
    nml_mat_free(LUP->L);
    nml_mat_free(LUP->U);
    nml_mat_free(LUP->P);
    free(LUP);
}

nml_mat_lup* nml_mat_LU(nml_mat *matrix){
	// Any square matrix can be written as LU
	// After decomposing A, easy to solve Ax = b
	// LUx = b, where Ux = y
	// Ly = b, which is solved using forward substitution
	// solving x via back substitution(Ux = y)
	

	//DOOLITES => L = identity with lower triangle


    nml_mat* copy = nml_mat_cp(matrix);
	unsigned int i,j,num_permutations;	
	int pivot;
	if(!matrix->isSquare){
		perror("Matrix is not square");
	}
	nml_mat *L = nml_mat_new(matrix->num_rows,matrix->num_rows);
	nml_mat *U = nml_mat_cp(matrix);
	nml_mat *P = nml_mat_iden(matrix->num_rows);
	for(j = 0;j<U->num_cols;++j){
		pivot = nml_mat_get_col_pivot(matrix,j,j);
	}

}


