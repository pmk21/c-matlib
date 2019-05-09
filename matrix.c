#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
    free(matrix);
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

Matrix* matrix_clone(Matrix *matrix)
{
    Matrix *mat;
    int rows = matrix->rows,
        cols = matrix->cols;

    /* Allocate space for the matrix structure */
    mat = (Matrix *)malloc(sizeof(Matrix));

    /* Initialize all the values */
    mat->rows = rows;
    mat->cols = cols;

    mat->data = (double **)malloc(sizeof(double *) * rows);

    for (int i=0;i<rows;i++)
    {
        mat->data[i] = (double *)malloc(sizeof(double) * cols);

        for (int j=0;j<cols;j++)
        {
            mat->data[i][j] = matrix->data[i][j];    
        }
    }

    return mat;
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
        printf("Warning: Orders of matrices do not match returning NULL!\n");
        return NULL;
    }
    else
    {
        /* Performs matrix addition */
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

Matrix* matrix_sub(Matrix *matrixA, Matrix *matrixB)
{
    /* Checks if the matrices are of the same order */
    if ( (matrixA->rows != matrixB->rows) || (matrixA->cols != matrixB->cols) ) 
    {
        printf("Warning: Orders of matrices do not match returning NULL!\n");
        return NULL;
    }
    else
    {   
        /* Performs matrixA - matrixB */
        Matrix *matrixC = (Matrix *)malloc(sizeof(Matrix));
        matrixC->rows = matrixA->rows;
        matrixC->cols = matrixA->cols;

        matrixC->data = (double **)malloc(sizeof(double *) * matrixC->rows);
        
        for (int i=0;i<matrixC->rows;i++)
        {
            matrixC->data[i] = (double *)malloc(sizeof(double) * matrixC->cols);
            
            for (int j=0;j<matrixC->cols;j++)
            {
                matrixC->data[i][j] = matrixA->data[i][j] - matrixB->data[i][j]; 
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

int matrix_det(Matrix *matrix)
{
    if (matrix->cols != matrix->rows)
    {    
        printf("Warning: Determinant exists only for square matrices!");
    }
    else
    {
        Matrix *mat;
        mat = matrix_clone(matrix);
        int swaps = 0;
        double det = 1.0, factor;

        for (int i=0;i<mat->rows;i++)
        {
            if (mat->data[i][i] == 0)
            {
                int k;
                for (k=i+1;k<mat->rows;k++)
                {
                    if (mat->data[k][i])
                        break;
                }

                if (k != mat->rows)
                {   
                    swap_row(mat, k, i);
                    swaps++;
                }
            }
            
            for (int j=i+1;j<mat->rows;j++)
            {
                if (mat->data[i][i] == 0)
                    continue;
                
                factor = mat->data[j][i]/mat->data[i][i];
                
                if (factor == 0)
                    continue;
                
                reduce(mat, i, j, factor);
            }

            matrix_print(mat);
        }

        for (int i=0;i<mat->rows;i++)
        {
            det = det * mat->data[i][i];
        }

        matrix_free(mat);

        return det * pow(-1, swaps);
    }
    
}

void reduce(Matrix* matrix, int i, int j, double factor)
{
    if (matrix->rows < i || matrix->rows < j)
        return;
    
    for (int k=0;k<matrix->cols;k++)
    {
        matrix->data[j][k] -= factor * matrix->data[i][k];
    }
}

void swap_row(Matrix* matrix, int i, int j)
{
    double temp = 0.0;

    for (int k=0; k<matrix->rows; k++) 
    { 
        temp = matrix->data[i][k];
        matrix->data[i][k] = matrix->data[j][k];
        matrix->data[j][k] = temp;
    } 
}