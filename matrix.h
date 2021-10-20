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
double nml_rand_interval(int min,int max);

//method to create a random matrix
nml_mat* nml_mat_rnd(unsigned int num_rows, unsigned int num_cols,int min,int max);

//method to create a square matrix
nml_mat* nml_mat_sqr(unsigned int rowCol);

//method to create an identity matrix
nml_mat* nml_mat_iden(unsigned int rowCol);

//method to copy matrix to another matrix
nml_mat* nml_mat_cp(nml_mat* matrix);

//read matrix from file[FIX FOR WINDOWS]
nml_mat* nml_mat_fromfile(FILE* f);
//FILE* pointer to the file which we are using


//method to create a random square matrix
nml_mat* nml_mat_sqr_rnd(unsigned int rowCol,int min,int max);

//check if matrix has same dimensions
int nml_mat_equaldim(nml_mat* m1,nml_mat* m2);

//check matrix equality
int nml_mat_equal(nml_mat* m1,nml_mat* m2);

#endif
