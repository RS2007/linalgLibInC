#ifndef MATRIX_H
#define MATRIX_H

typedef struct nml_mat_s
{
	unsigned int num_rows;
	unsigned int num_cols;
	double **data;
	int isSquare;
} nml_mat;

//constructing a matrix frame(allocating memory -> all elements are 0)
nml_mat *nml_mat_new(unsigned int num_rows, unsigned int num_cols);

//destructor (destroys memory -> each malloc should have its own free)
void nml_mat_free(nml_mat *matrix);

//generate random number in interval
double nml_rand_interval(int min, int max);

//method to create a random matrix
nml_mat *nml_mat_rnd(unsigned int num_rows, unsigned int num_cols, int min, int max);

//method to create a square matrix
nml_mat *nml_mat_sqr(unsigned int rowCol);

//method to create an identity matrix
nml_mat *nml_mat_iden(unsigned int rowCol);

//method to copy matrix to another matrix
nml_mat *nml_mat_cp(nml_mat *matrix);

//read matrix from file[FIX FOR WINDOWS]
nml_mat *nml_mat_fromfile(FILE *f);
//FILE* pointer to the file which we are using

//method to create a random square matrix
nml_mat *nml_mat_sqr_rnd(unsigned int rowCol, int min, int max);

//check if matrix has same dimensions
int nml_mat_equaldim(nml_mat *m1, nml_mat *m2);

//check matrix equality
int nml_mat_equal(nml_mat *m1, nml_mat *m2);

//retreiving the column from the matrixx
nml_mat *nml_mat_get_column(nml_mat *matrix, unsigned int col);

//retreiving the row from  the matrixx
nml_mat *nml_mat_get_row(nml_mat *matrix, unsigned int row);

//setting all elements in matrix to a single number
nml_mat *nml_set_all_elements(nml_mat *matrix, double num);

//setting all main diagonal elements to a single number
nml_mat *nml_set_diagonal_elements(nml_mat *matrix, double num);

//multiply a row with a scalar
nml_mat *nml_row_multipy_scalar(nml_mat *matrix, unsigned int row, int scalar);

//multiply a col with a scalar
nml_mat *nml_col_multiply_scalar(nml_mat *matrix, unsigned int col, int scalar);

//adding rows in a matrix(for gaussian elimination or something like that
nml_mat *nml_rows_add(nml_mat *matrix, unsigned int addendum, unsigned int original, int multiplier);

//printing a matrix
void nml_mat_print(nml_mat *matrix);

//multipy matrix with a scalar
nml_mat *nml_mat_multiply_scalar(nml_mat *matrix, int scalar);

//removing a column in a matrix
nml_mat *nml_mat_remove_column(nml_mat *matrix, unsigned int column);

//removing a row in a matrix
nml_mat *nml_mat_remove_row(nml_mat *matrix, unsigned int row);

//swapping rows
nml_mat* nml_mat_swap_row(nml_mat* matrix,unsigned int row1,unsigned int row2);

//swapping columns
nml_mat* nml_mat_swap_column(nml_mat* matrix,unsigned int col1,unsigned int col2);

//horizontal concatenation of two matrices
nml_mat* nml_mat_horizontal_concat(nml_mat* matrix1,nml_mat* matrix2);

//vertical concatenation
nml_mat* nml_mat_vertical_concat(nml_mat* matrix1,nml_mat* matrix2);

//add matrices
nml_mat* nml_mat_add(nml_mat* matrix1,nml_mat* matrix2);

//subtracting matrices

//multiplying matrices

//row echelon form

//reduced row echelon form

// LU(P) decomposition

//forward substituition linear system

//backward substitution linear system

//linear system using LU(P)

//inverse of matrix using LU(P)

//determinant of matrix using LU(P)

//QR decomposition

#endif
