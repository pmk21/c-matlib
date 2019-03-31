// Initialize functions for operations on matrices

/* Matrix structure */
struct matrix{
    int rows;   // Number of rows
    int cols;   // Number of columns
    double **data; // Data of a matrix
} matrix;

typedef struct matrix Matrix;

/* Initialize integer matrix */
Matrix* matrix_init();

/* Get value at specified index */
int matrix_index(Matrix *matrix, int row, int col);

/* Free matrix */
void matrix_free(Matrix *matrix);

/* Print the matrix */
void matrix_print(Matrix *matrix);

// Basic operations

/* Transposes the matrix */
void matrix_T(Matrix *matrix);

/* Matrix Addition */
Matrix* matrix_add(Matrix *matrixA, Matrix *matrixB);

/* Matrix Subtraction */
Matrix* matrix_sub(Matrix *matrixA, Matrix *matrixB);

/* Scalar multiplication */
void matrix_smul(Matrix *matrix, double scalar);