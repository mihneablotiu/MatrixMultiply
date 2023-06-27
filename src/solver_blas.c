/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"
#include "cblas.h"
#include <string.h>

/* 
 * Add your BLAS implementation here
 */
double* my_solver(int N, double *A, double *B) {
	int number_of_elements = N * N;

	double *result_matrix = malloc(number_of_elements * sizeof(*result_matrix));
	if (result_matrix == NULL) {
		fprintf(stderr, "%s\n", "Error at malloc for the partial results matrix");
		exit(1);
	}

	/* Doing the B^T * B^T operation firstly and the result being saved in the
	result_matrix */
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
	 			N, N, N, 1.0, B, N, B, N, 0.0, result_matrix, N);

	/* Doing the A * B operation taking into consideration the fact that A is a
	upper triangular matrix and saving the result in B*/	
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans,
			    CblasNonUnit, N, N, 1.0, A, N, B, N);

	/* Doing the multiplication between the last result saved in B and A^T taking
	into consideration that A is a upper triangular matrix so that means that A^T
	is a lower triangular and saving the result in B */
	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasTrans,
				CblasNonUnit, N, N, 1.0, A, N, B, N);

	/* Doing the add between the previous two results (stored in B and in result_matrix)
	and save the result in the result_matrix which we return from the function */
	cblas_daxpy(number_of_elements, 1.0, B, 1, result_matrix, 1);

	return result_matrix;
}
