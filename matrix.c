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

Matrix* matrix_create(int row, int col)
{
    Matrix *matrix;
    
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
            matrix->data[i][j] = 0;
        }
    }
    return matrix;
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
Matrix* matrix_mul(Matrix *matrixA, Matrix *matrixB)
{
    if (matrixA->cols != matrixB->rows)
    {
        printf("Warning: Multiplication not possible on matrices of the given order, returning NULL!\n");
        return NULL;
    }
    else if((matrixA->cols==matrixB->cols && matrixA->rows==matrixB->rows) && (matrixA->rows)%2 == 0 && (matrixB->rows)%2 == 0 && (matrixB->cols)%2 == 0 && (matrixB->cols)%2 == 0)
    {
        /*Performs Strassens's Matrix Multiplication*/
        Matrix *matrixC = (Matrix *)malloc(sizeof(Matrix));
        matrixC->rows = matrixA->rows;
        matrixC->cols = matrixA->cols;
        int n = matrixA->rows;
        matrixC->data = (double **)malloc(sizeof(double *) * matrixC->rows);
        for (int i=0;i<matrixC->rows;i++)
        {
            matrixC->data[i] = (double *)malloc(sizeof(double) * matrixC->cols);
            
            for (int j=0;j<matrixC->cols;j++)
            {
                matrixC->data[i][j] = 0; 
            }
        }
        matrixC = matrix_mul_strassens(matrixA, matrixB, n);
        return matrixC;
    }
    else{
        Matrix *matrixC = (Matrix *)malloc(sizeof(Matrix));
        matrixC->rows = matrixA->rows;
        matrixC->cols = matrixB->cols;
        matrixC->data = (double **)malloc(sizeof(double *) * matrixC->rows);
        for (int i=0;i<matrixC->rows;i++)
        {
            matrixC->data[i] = (double *)malloc(sizeof(double) * matrixC->cols);
            
            for (int j=0;j<matrixC->cols;j++)
            {
                matrixC->data[i][j] = 0; 
            }
        }
        int r1 = matrixA->rows;
        int c1 = matrixB->rows;
        int c2 = matrixB->cols;
        for(int i=0; i<r1; ++i){
            for(int j=0; j<c2; ++j){
                for(int k=0; k<c1; ++k){
                    matrixC->data[i][j] += (matrixA->data[i][k])*(matrixB->data[k][j]);
                }
            }
        }
        return matrixC;
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
        //TODO: Test this function more
        
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

void matrix_LU(Matrix *matrix, Matrix **l, Matrix **u)
{
    // TODO: Add conditions for checking whether LU decomposition exists
    // TODO: Implement LU decomposition for rectangular matrices

    if (matrix->rows != matrix->cols)
    {
        printf("ERROR: LU Decomposition for rectangular matrices not implemented\n");
        return;
    }

    Matrix *L = matrix_create(matrix->rows, matrix->rows);
    Matrix *U = matrix_create(matrix->rows, matrix->rows);

    for (int i = 0; i < matrix->rows; i++)
    {
        // Upper Triangular
        for (int k = i; k < matrix->rows; k++)
        { 
            // Summation of L(i, j) * U(j, k) 
            int sum = 0; 
            for (int j = 0; j < i; j++) 
                sum += (L->data[i][j] * U->data[j][k]); 
  
            // Evaluating U(i, k) 
            U->data[i][k] = matrix->data[i][k] - sum; 
        }
  
        // Lower Triangular 
        for (int k = i; k < matrix->rows; k++) 
        {
            
            if (i == k)
                L->data[i][i] = 1; // Diagonal as 1
            else 
            {
                // Summation of L(k, j) * U(j, i)
                int sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (L->data[k][j] * U->data[j][i]);
  
                // Evaluating L(k, i)
                L->data[k][i] = (matrix->data[k][i] - sum) / U->data[i][i];
            }
        }
    }

    *l = L;
    *u = U;
}
Matrix* matrix_mul_strassens(Matrix *matrixA, Matrix *matrixB, int n){

    if(n<=2){
        /* Multiply the matrices if they are of the order 2X2 */
        Matrix *matrixD = (Matrix *)malloc(sizeof(Matrix));
        matrixD->rows = 2;
        matrixD->cols = 2;
        int n = matrixA->rows;
        matrixD->data = (double **)malloc(sizeof(double *) * matrixD->rows);
        for (int i=0;i<matrixD->rows;i++)
        {
            matrixD->data[i] = (double *)malloc(sizeof(double) * matrixD->cols);     
            for (int j=0;j<matrixD->cols;j++)
            {
                matrixD->data[i][j] = 0; 
            }
        }
        double m1, m2, m3, m4 , m5, m6, m7;
        //m1= (a[0][0] + a[1][1]) * (b[0][0] + b[1][1]);
        m1 = (matrixA->data[0][0] + matrixA->data[1][1]) * (matrixB->data[0][0] + matrixB->data[1][1]);
        //m2= (a[1][0] + a[1][1]) * b[0][0];
        m2 = (matrixA->data[1][0] + matrixA->data[1][1]) * matrixB->data[0][0];
        //m3= a[0][0] * (b[0][1] - b[1][1]);
        m3 = matrixA->data[0][0] * (matrixB->data[0][1] - matrixB->data[1][1]);
        //m4= a[1][1] * (b[1][0] - b[0][0]);
        m4 = matrixA->data[1][1] * (matrixB->data[1][0] - matrixB->data[0][0]);
        //m5= (a[0][0] + a[0][1]) * b[1][1];
        m5 = (matrixA->data[0][0] + matrixA->data[0][1]) * matrixB->data[1][1];
        //m6= (a[1][0] - a[0][0]) * (b[0][0]+b[0][1]);
        m6 = (matrixA->data[1][0] - matrixA->data[0][0]) * (matrixB->data[0][0] + matrixB->data[0][1]);
        //m7= (a[0][1] - a[1][1]) * (b[1][0]+b[1][1]);
        m7 = (matrixA->data[0][1] - matrixA->data[1][1]) * (matrixB->data[1][0] + matrixB->data[1][1]);

        //c[0][0] = m1 + m4- m5 + m7;
        matrixD->data[0][0] = m1 + m4 - m5 + m7;
        //c[0][1] = m3 + m5;
        matrixD->data[0][1] = m3 + m5;
        //c[1][0] = m2 + m4;
        matrixD->data[1][0] = m2 + m4;
        //c[1][1] = m1 - m2 + m3 + m6;
        matrixD->data[1][1] = m1 - m2 + m3 + m6;
        return matrixD;
    }
    else{
        int index1=0;
        int index2=0;
        int n = matrixA->rows;
        Matrix *matrixA11 = (Matrix*)malloc(sizeof(Matrix));
        Matrix *matrixA12 = (Matrix*)malloc(sizeof(Matrix));
        Matrix *matrixA21 = (Matrix*)malloc(sizeof(Matrix));
        Matrix *matrixA22 = (Matrix*)malloc(sizeof(Matrix));

        Matrix *matrixB11 = (Matrix*)malloc(sizeof(Matrix));
        Matrix *matrixB12 = (Matrix*)malloc(sizeof(Matrix));
        Matrix *matrixB21 = (Matrix*)malloc(sizeof(Matrix));
        Matrix *matrixB22 = (Matrix*)malloc(sizeof(Matrix));

        matrixA11->rows = n/2;
        matrixA11->cols = n/2;
        matrixA12->rows = n/2;
        matrixA12->cols = n/2;
        matrixA21->rows = n/2;
        matrixA21->cols = n/2;
        matrixA22->rows = n/2;
        matrixA22->cols = n/2;

        matrixB11->rows = n/2;
        matrixB11->cols = n/2;
        matrixB12->rows = n/2;
        matrixB12->cols = n/2;
        matrixB21->rows = n/2;
        matrixB21->cols = n/2;
        matrixB22->rows = n/2;
        matrixB22->cols = n/2;

        matrixA11->data = (double **)malloc(sizeof(double *) * matrixA->rows/2);
        matrixA12->data = (double **)malloc(sizeof(double *) * matrixA->rows/2);
        matrixA21->data = (double **)malloc(sizeof(double *) * matrixA->rows/2);
        matrixA22->data = (double **)malloc(sizeof(double *) * matrixA->rows/2);

        matrixB11->data = (double **)malloc(sizeof(double *) * matrixA->rows/2);
        matrixB12->data = (double **)malloc(sizeof(double *) * matrixA->rows/2);
        matrixB21->data = (double **)malloc(sizeof(double *) * matrixA->rows/2);
        matrixB22->data = (double **)malloc(sizeof(double *) * matrixA->rows/2);
        
        for(int j=0; j<n/2; j++){
            matrixA11->data[j] = (double *)malloc(sizeof(double)*n/2);
            matrixB11->data[j] = (double *)malloc(sizeof(double)*n/2);
            for(int k=0; k<n/2; k++){
                matrixA11->data[j][k] = matrixA->data[index1][index2];
                matrixB11->data[j][k] = matrixB->data[index1][index2];
                index2++;
            }
            index1++;
            index2=0;
        }
        index1 = n/2;
        index2 = 0;
        for(int j=0; j<n/2; j++){
            matrixA21->data[j] = (double *)malloc(sizeof(double)*n/2);
            matrixB21->data[j] = (double *)malloc(sizeof(double)*n/2);
            for(int k=0; k<n/2; k++){
                matrixA21->data[j][k] = matrixA->data[index1][index2];
                matrixB21->data[j][k] = matrixB->data[index1][index2];
                index2++;
            }
            index1++;
            index2=0;
        }
        index1 = 0;
        index2 = n/2;
        for(int j=0; j<n/2; j++){
            matrixA12->data[j] = (double *)malloc(sizeof(double)*n/2);
            matrixB12->data[j] = (double *)malloc(sizeof(double)*n/2);
            for(int k=0; k<n/2; k++){
                matrixA12->data[j][k] = matrixA->data[index1][index2];
                matrixB12->data[j][k] = matrixB->data[index1][index2];
                index2++;
            }
            index1++;
            index2=n/2;
        }
        index1 = n/2;
        index2 = n/2;
        for(int j=0; j<n/2; j++){
            matrixA22->data[j] = (double *)malloc(sizeof(double)*n/2);
            matrixB22->data[j] = (double *)malloc(sizeof(double)*n/2);
            for(int k=0; k<n/2; k++){
                matrixA22->data[j][k] = matrixA->data[index1][index2];
                matrixB22->data[j][k] = matrixB->data[index1][index2];
                index2++;
            }
            index1++;
            index2=n/2;
        }

        /*Making the Recursive Calls*/
        Matrix *mat1;
        Matrix *mat2;
        
        Matrix *sum1;
        mat1 = matrix_mul_strassens(matrixA11,matrixB11,n/2);
        mat2 = matrix_mul_strassens(matrixA12,matrixB21,n/2);
        sum1 = matrix_add(mat1,mat2);
        
        Matrix *sum2;
        mat1 = matrix_mul_strassens(matrixA11,matrixB12,n/2);
        mat2 = matrix_mul_strassens(matrixA12,matrixB22,n/2);
        sum2 = matrix_add(mat1,mat2);
        
        Matrix *sum3;
        mat1 = matrix_mul_strassens(matrixA21,matrixB11,n/2);
        mat2 = matrix_mul_strassens(matrixA22,matrixB21,n/2);
        sum3 = matrix_add(mat1,mat2);
        
        Matrix *sum4;
        mat1 = matrix_mul_strassens(matrixA21,matrixB12,n/2);
        mat2 = matrix_mul_strassens(matrixA22,matrixB22,n/2);
        sum4 = matrix_add(mat1,mat2);
        
        Matrix *matrixE = (Matrix *)malloc(sizeof(Matrix));
        matrixE->rows = n;
        matrixE->cols = n;
        matrixE->data = (double **)malloc(sizeof(double *) * matrixE->rows);
        for (int i=0;i<matrixE->rows;i++)
        {
            matrixE->data[i] = (double *)malloc(sizeof(double) * matrixE->cols);
            
            for (int j=0;j<matrixE->cols;j++)
            {
                matrixE->data[i][j] = 0; 
            }
        }
        int var1=0;
        int var2=0;
        for (int i=0;i<sum1->rows;i++)
        {
            for (int j=0;j<sum1->cols;j++)
            {
                matrixE->data[var1][var2] = sum1->data[i][j]; 
                var2++;
            }
            var2=0;
            var1++;
        }

        var1 = 0;
        var2 = n/2;
        for (int i=0;i<sum2->rows;i++)
        {
            for (int j=0;j<sum2->cols;j++)
            {
                matrixE->data[var1][var2] = sum2->data[i][j]; 
                var2++;
            }
            var2=n/2;
            var1++;
        }
        
        var1 = n/2;
        var2 = 0;
        for (int i=0;i<sum3->rows;i++)
        {
            for (int j=0;j<sum3->cols;j++)
            {
                matrixE->data[var1][var2] = sum3->data[i][j]; 
                var2++;
            }
            var2=0;
            var1++;
        }
        var1 = n/2;
        var2 = n/2;
        for (int i=0;i<sum4->rows;i++)
        {
            for (int j=0;j<sum4->cols;j++)
            {
                matrixE->data[var1][var2] = sum4->data[i][j]; 
                var2++;
            }
            var2=n/2;
            var1++;
        }
        return matrixE;
    }
}