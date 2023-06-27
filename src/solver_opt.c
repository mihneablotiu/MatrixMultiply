/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

/* The function that multiplies a upper triangular matrix to a random matrix from
the left, taking into consideration that the left matrix has this special property
so we do not compute unnecessary operations using pointers to access vectors,
register variables for additional optimizations and reordering the loops in
i-k-j order to have a better memory access  */
void multiply_upper_triangular_with_random(double *upper_triangular, double *random, double *result, int N)
{
	int i, j, k;
	register double *orig_upper_triangular = upper_triangular;
	register double *orig_random = random;
	register double *orig_result = result;

	for (i = 0; i < N; i++) {
		orig_random = random;

		for (k = i, orig_upper_triangular += k, orig_random += k * N; k < N; k++) {
			orig_result = result + N * i;

			for (j = 0; j < N; j++) {
				*orig_result += (*orig_upper_triangular) * (*orig_random);

				orig_random++;
				orig_result++;
			}

			orig_upper_triangular++;
		}
	}
}

/* The function that multiplies a lower triangular matrix to a random matrix from
the right, taking into consideration that the right matrix has this special property
so we do not compute unnecessary operations using pointers to access vectors and
register variables for additional optimizations */
void multiply_random_with_lower_triangular(double *random, double *lower_triangular, double *result, int N)
{
	int i, j, k;
	double *orig_lower_triangular;
	register double *current_random, *orig_random;
	register double current_sum;

	for (i = 0; i < N; i++) {
		orig_random = random + N * i;

		for (j = 0; j < N; j++) {
			current_random = orig_random;
			orig_lower_triangular = lower_triangular + j;

			current_sum = 0.0;

			for (k = j, current_random += k, orig_lower_triangular += k * N; k < N; k++) {
				current_sum += (*current_random) * (*orig_lower_triangular);

				current_random++;
				orig_lower_triangular += N;
			}

			*(result + i * N + j) = current_sum;
		}
	}
}

/* The function that clasically multiplies i-j-k two random matrixes without any
particular properties using pointers to access vectors, register variables
for additional optimizations and reordering the loops in i-k-j order to have a
better memory access */
void multiply_random_with_random(double *A, double *B, double *result, int N)
{
	int i, j, k;
	register double *orig_A = A;
	register double *orig_B = B;
	register double *orig_result = result;

	for (i = 0; i < N; i++) {
		orig_B = B;

		for (k = 0; k < N; k++) {
			orig_result = result + N * i;

			for (j = 0; j < N; j++) {
				*orig_result += (*orig_A) * (*orig_B);

				orig_B++;
				orig_result++;
			}

			orig_A++;
		}
	}
}

/* The function that just computes the sum of two random matrixes using pointers to access
vectors and register variables for additional optimizations */
void add_random_with_random(double *A, double *B, double *result, int N)
{
	register double *orig_A = A;
	register double *orig_B = B;
	register double *orig_result = result;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			*orig_result = *orig_A + *orig_B;

			orig_result++;
			orig_A++;
			orig_B++;
		}
	}
}

/* The function that computes the transpose of a matrix in-place (without additional)
memory by switching the i-j element with the j-i element using pointers to access
vectors and register variables for additional optimizations */
void transpose_matrix(double *matrix, int N)
{
	register double *orig_A;
	register double *orig_B;

	int i = 0, j = 0;

	for (i = 0; i < N; i++) {
		orig_A = matrix + N * i;
		orig_B = matrix + i;

		for (j = 0; j < i; j++) {
			double aux = *orig_A;
			*orig_A = *orig_B;
			*orig_B = aux;

			orig_B += N;
			orig_A++;
		}
	}

}

/* The function that checks whether a malloc failed or not */
void malloc_check(double *matrix)
{
	if (matrix == NULL) {
		fprintf(stderr, "%s\n", "Error at malloc for the partial results matrix");
		exit(1);
	}
}

void matrix_alloc(double **AB, double **ABA_transpose,
				  double **B_transposeB_transpose, double **final_result,
				  int number_of_elements)
{
	*AB = calloc(number_of_elements, sizeof(**AB));
	malloc_check(*AB);

	*ABA_transpose = calloc(number_of_elements, sizeof(**ABA_transpose));
	malloc_check(*ABA_transpose);

	*B_transposeB_transpose = calloc(number_of_elements, sizeof(**B_transposeB_transpose));
	malloc_check(*B_transposeB_transpose);

	*final_result = malloc(number_of_elements * sizeof(**final_result));
	malloc_check(*final_result);
}

void free_memory(double **AB, double **ABA_transpose, double **B_transposeB_transpose)
{
	free(*AB);
	free(*ABA_transpose);
	free(*B_transposeB_transpose);

	AB = NULL;
	ABA_transpose = NULL;
	B_transposeB_transpose = NULL;
}

/*
 * Add your optimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	int number_of_elements = N * N;

	double *AB, *ABA_transpose;
	double *B_transposeB_transpose;
	double *final_result;

	/* We reserve memory for each of the intermediate steps of the operation*/
	matrix_alloc(&AB, &ABA_transpose, &B_transposeB_transpose,
				 &final_result, number_of_elements);

	/* We multiply A with B and save the result in AB */
	multiply_upper_triangular_with_random(A, B, AB, N);

	/* We transpose the matrix A and now A is A^T */
	transpose_matrix(A, N);

	/* We multiply the AB matrix stored before with A^T and store the result
	in ABA_Transpose */
	multiply_random_with_lower_triangular(AB, A, ABA_transpose, N);

	/* We transpose the matrix B and now B is B^T */
	transpose_matrix(B, N);

	/* We multiply B with B which now is actually B^T * B^T and save the result
	in B_transposeB_transpose */
	multiply_random_with_random(B, B, B_transposeB_transpose, N);

	/* We add the previous two results and store the final result*/
	add_random_with_random(ABA_transpose, B_transposeB_transpose, final_result, N);

	/* We free the memory of each intermediate matrix and return the final result */
	free_memory(&AB, &ABA_transpose, &B_transposeB_transpose);
	return final_result;
}