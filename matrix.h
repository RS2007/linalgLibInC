#ifndef MATRIX_H
#define MATRIX_H

typedef struct nml_mat_s
{
	unsigned int num_rows;
	unsigned int num_cols;
	double **data;
	int isSquare;
} nml_mat;

typedef struct nml_mat_lup{
    nml_mat* L;
    nml_mat* U;
    nml_mat* P;
    unsigned int num_permutations;
} nml_mat_lup;


#define NML_MIN_COEFF 0.00001

// constructing a matrix frame(allocating memory -> all elements are 0)
nml_mat *nml_mat_new(unsigned int num_rows, unsigned int num_cols);

// destructor (destroys memory -> each malloc should have its own free)
void nml_mat_free(nml_mat *matrix);

// constructing the LUP frame(allocating memory and all the elements are 0)
nml_mat_lup *nml_mat_lup_new(nml_mat* L,nml_mat* U, nml_mat* P,unsigned int num_permutations);

// destructor (destroys memory->each malloc should have its own free(RULE))
void nml_mat_lup_free(nml_mat_lup *LUP);

// generate random number in interval
double nml_rand_interval(int min, int max);

// method to create a random matrix
nml_mat *nml_mat_rnd(unsigned int num_rows, unsigned int num_cols, int min, int max);

// method to create a square matrix
nml_mat *nml_mat_sqr(unsigned int rowCol);

// method to create an identity matrix
nml_mat *nml_mat_iden(unsigned int rowCol);

// method to copy matrix to another matrix
nml_mat *nml_mat_cp(nml_mat *matrix);

// read matrix from file[FIX FOR WINDOWS]
nml_mat *nml_mat_fromfile(FILE *f);
// FILE* pointer to the file which we are using

// method to create a random square matrix
nml_mat *nml_mat_sqr_rnd(unsigned int rowCol, int min, int max);

// check if matrix has same dimensions
int nml_mat_equaldim(nml_mat *m1, nml_mat *m2);

// check matrix equality
int nml_mat_equal(nml_mat *m1, nml_mat *m2);

// retreiving the column from the matrixx
nml_mat *nml_mat_get_column(nml_mat *matrix, unsigned int col);

// retreiving the row from  the matrixx
nml_mat *nml_mat_get_row(nml_mat *matrix, unsigned int row);

// setting all elements in matrix to a single number
nml_mat *nml_set_all_elements(nml_mat *matrix, double num);

// setting all main diagonal elements to a single number
nml_mat *nml_set_diagonal_elements(nml_mat *matrix, double num);

// multiply a row with a scalar
nml_mat *nml_row_multipy_scalar(nml_mat *matrix, unsigned int row, double scalar);

// multiply a col with a scalar
nml_mat *nml_col_multiply_scalar(nml_mat *matrix, unsigned int col, double scalar);

// adding rows in a matrix(for gaussian elimination or something like that
nml_mat *nml_rows_add(nml_mat *matrix, unsigned int addendum, unsigned int original, double multiplier);

// printing a matrix
void nml_mat_print(nml_mat *matrix);

// multipy matrix with a scalar
nml_mat *nml_mat_multiply_scalar(nml_mat *matrix, int scalar);

// removing a column in a matrix
nml_mat *nml_mat_remove_column(nml_mat *matrix, unsigned int column);

// removing a row in a matrix
nml_mat *nml_mat_remove_row(nml_mat *matrix, unsigned int row);

// swapping rows
nml_mat *nml_mat_swap_row(nml_mat *matrix, unsigned int row1, unsigned int row2);

// swapping columns
nml_mat *nml_mat_swap_column(nml_mat *matrix, unsigned int col1, unsigned int col2);

// horizontal concatenation of two matrices
nml_mat *nml_mat_horizontal_concat(int mnum, nml_mat **marr);

// vertical concatenation
nml_mat *nml_mat_vertical_concat(int mnum, nml_mat **marr);

// add matrices
nml_mat *nml_mat_add(nml_mat *matrix1, nml_mat *matrix2);

// subtracting matrices
nml_mat *nml_mat_subtract(nml_mat *matrix1, nml_mat *matrix2);

// multiplying matrices(using O(n^3))
nml_mat *nml_mat_mul_naive(nml_mat *matrix1, nml_mat *matrix2);

// multiplyting matrices(using O(n^2.87))
nml_mat *nml_mat_mul_strassen(nml_mat *matrix1, nml_mat *matrix2);

// get the pivot in a column
int nml_mat_get_col_pivot(nml_mat *matrix, unsigned introw1, unsigned int row2);

// row echelon form
nml_mat *nml_mat_ref(nml_mat *matrix);

// reduced row echelon form
nml_mat *nml_mat_rref(nml_mat *matrix);

nml_mat* nml_mat_swap_row_LU(nml_mat* matrix,unsigned int row1,unsigned int row2);



// LU(P) decomposition
nml_mat_lup *nml_mat_LU(nml_mat *matrix);

// forward substituition linear system
nml_mat *solve_linear_forward();

// backward substitution linear system
nml_mat *solve_linear_backward();

// linear system using LU(P)
nml_mat *solve_linear_LU();

// inverse of matrix using LU(P)
nml_mat *nml_mat_inverse();

// determinant of matrix using LU(P)
nml_mat *nml_mat_determinant();

// QR decomposition
nml_mat *nml_mat_QR();

#endif
