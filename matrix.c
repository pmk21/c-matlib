#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

Matrix* matrix_init()
{
    Matrix *matrix;
    int row, col;
    double temp;

    scanf("%d", &row);
    scanf("%d", &col);

    /* Allocate space for the matrix structure */
    matrix = (Matrix *)malloc(sizeof(Matrix));

    /* Initialize all the values */
    matrix->rows = row;
    matrix->cols = col;

    matrix->data = (double **)malloc(sizeof(double *) * row);

    for (int i=0;i<row;i++)
    {
        matrix->data[i] = (double *)malloc(sizeof(double) * col);

        for (int j=0;j<col;j++)
        {
            scanf("%lf", &temp);
            matrix->data[i][j] = temp;    
        }
    }

    return matrix;
}


int matrix_index(Matrix *matrix, int row, int col)
{
    /* Get an element at a particular index */ 
    return matrix->data[row][col];
}


void matrix_free(Matrix *matrix)
{
    /* De-allocates memory for the main array holding the matrix */
    for (int i=0;i<matrix->rows;i++)
    {
        free(matrix->data[i]);
    }
    free(matrix->data);
}

void matrix_print(Matrix *matrix)
{

    if (matrix == NULL)
    {
        printf("Empty Matrix\n");
        return;
    }

    /* Print the order of the matrix */
    printf("Size: (%d X %d)\n", matrix->rows, matrix->cols);

    /* Prints the rows and columns */
    for (int i=0;i<matrix->rows;i++)
    {
        for (int j=0;j<matrix->cols;j++)
        {
            printf("%.2lf ", matrix->data[i][j]);
        }
        printf("\n");
    }
}

void matrix_T(Matrix *matrix)
{

    /* Check if it is a square matrix or not */
    if (matrix->rows == matrix->cols)
    {
        /* Straightforward for finding transpose of a square matrix */
        int temp;

        for (int i=0;i<(matrix->rows - 1);i++)
        {
            for (int j=i;j<(matrix->cols);j++)
            {
                temp = matrix->data[i][j];
                matrix->data[i][j] = matrix->data[j][i];
                matrix->data[j][i] = temp;
            }
        }

    }
    else
    {
        /* Finding transpose of a rectangular matrix */
        double **matT, temp;
        
        matT = (double **)malloc(sizeof(double *) * matrix->cols);
        
        for (int i=0;i<matrix->cols;i++)
        {    
            matT[i] = (double *)malloc(sizeof(double) * matrix->rows);
        }
        
        for (int i=0;i<matrix->rows;i++)
        {
            for (int j=0;j<matrix->cols;j++)
            {
                matT[j][i] = matrix->data[i][j]; 
            }
        }

        matrix_free(matrix);

        matrix->data = matT;
        temp = matrix->cols;
        matrix->cols = matrix->rows;
        matrix->rows = temp;
    }
}

Matrix* matrix_add(Matrix *matrixA, Matrix *matrixB)
{
    /* Checks if the matrices are of the same order */
    if ( (matrixA->rows != matrixB->rows) || (matrixA->cols != matrixB->cols) ) 
    {
        printf("Warning: Orders of matrices do not match\n");
        return NULL;
    }
    else
    {
        Matrix *matrixC = (Matrix *)malloc(sizeof(Matrix));
        matrixC->rows = matrixA->rows;
        matrixC->cols = matrixA->cols;

        matrixC->data = (double **)malloc(sizeof(double *) * matrixC->rows);
        
        for (int i=0;i<matrixC->rows;i++)
        {
            matrixC->data[i] = (double *)malloc(sizeof(double) * matrixC->cols);
            
            for (int j=0;j<matrixC->cols;j++)
            {
                matrixC->data[i][j] = matrixA->data[i][j] + matrixB->data[i][j]; 
            }
        }

        return matrixC;
    }
}

void matrix_smul(Matrix *matrix, double scalar)
{
    /* Performs scalar multiplication on a matrix */
    for (int i=0;i<matrix->rows;i++)
    {
        for (int j=0;j<matrix->cols;j++)
        {
            matrix->data[i][j] = scalar * matrix->data[i][j];
        }
    }
}