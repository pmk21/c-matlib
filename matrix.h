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

/* Create a copy of the matrix */
Matrix* matrix_clone(Matrix *matrix);

/* Return a Matrix structure given row and column size filled with zeros */
Matrix* matrix_create(int row, int col);

// Basic operations

/* Transposes the matrix */
void matrix_T(Matrix *matrix);

/* Matrix Addition */
Matrix* matrix_add(Matrix *matrixA, Matrix *matrixB);

/* Matrix Subtraction */
Matrix* matrix_sub(Matrix *matrixA, Matrix *matrixB);

/* Scalar Multiplication */
void matrix_smul(Matrix *matrix, double scalar);

/* Strassen's Multiplication */
Matrix* matrix_mul_strassens(Matrix *matrixA, Matrix *matrixB, int n);

/*Matrix Multiplication*/
Matrix* matrix_mul(Matrix *matrixA, Matrix *matrixB);


/* Scalar multiplication */
void matrix_smul(Matrix *matrix, double scalar);

/* Find determinant */
int matrix_det(Matrix *matrix);

/* Swap rows numbered i and j */
void swap_row(Matrix *matrix, int i, int j);

/* Reduce row i by factor times row j */
void reduce(Matrix *matrix, int i, int j, double factor);

/* 
    Perform LU decomposition 
    
    IMPORTANT: L, U inputs should be pointer-to-pointer-to-structure
    (i.e double pointer)

*/
void matrix_LU(Matrix *matrix, Matrix **L, Matrix **U);
