/* File     : matrix_utils.cpp
 * Author   : Michael Hicks
 * 
 * Description:
 * Contains utility functions for matrix operations. Primarily for
 * use in optimizing ECA_CD cancellation.
 */


#include "matrix_utils.h"

#define real_prod(a, b, c, d) (a*c - b*d)
#define comp_prod(a, b, c, d) (a*d + b*c)

/*  CofactorMatrix_3x3:
 *  Returns the matrix of cofactors of the 3x3 COMPLEX matrix A
 *  A[0] = 3x3 real data 
 *  A[1] = 3x3 imag data
 *  i.e. returns C such that adj(A) = C^(T)
 */
void CofactorMatrix_3x3(double A[][3][3], double C[][3][3]){
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            // Real component
            C[0][i][j] = (
                      (A[0][(i+1)%3][(j+1)%3] * A[0][(i+2)%3][(j+2)%3] 
                    -  A[1][(i+1)%3][(j+1)%3] * A[1][(i+2)%3][(j+2)%3]) 
                    - (A[0][(i+1)%3][(j+2)%3] * A[0][(i+2)%3][(j+1)%3] 
                    -  A[1][(i+1)%3][(j+2)%3] * A[1][(i+2)%3][(j+1)%3])
                    );
            
            // Complex Component
            C[1][i][j] = (
                      (A[0][(i+1)%3][(j+1)%3]*A[1][(i+2)%3][(j+2)%3] 
                    +  A[1][(i+1)%3][(j+1)%3]*A[0][(i+2)%3][(j+2)%3]) 
                    - (A[0][(i+1)%3][(j+2)%3]*A[1][(i+2)%3][(j+1)%3] 
                    +  A[1][(i+1)%3][(j+2)%3]*A[0][(i+2)%3][(j+1)%3])
                    );
        }
    }
}

/* Determinant_3x3:
 * Returns the determinant of the 3x3 COMPLEX matrix A
 * C = pre-computed matrix of cofactors
 * A[0], C[0] = 3x3 real data
 * A[1], C[1] = 3x3 imag data
 * det[0] = real component of determinant 
 * det[1] = imag component of determinant
 */
void Determinant_3x3(double A[][3][3], double C[][3][3], double det[2]){
    det[0] = 0;
    det[1] = 1;
         
    for(int i = 0; i < 3; i++){
        det[0] = det[0] + real_prod(A[0][0][i], A[1][0][i], 
                                    C[0][0][i], C[1][0][i]);
        
        det[1] = det[1] + comp_prod(A[0][0][i], A[1][0][i], 
                                    C[0][0][i], C[1][0][i]);
    }
}

/* InvertMatrix_3x3:
 * Returns the inverse of the 3x3 COMPLEX matrix A
 * A[0], B[0] = real data
 * A[1], B[1] = imag data
 * B = A^(-1)
 */
void InvertMatrix_3x3(double A[][3][3], double B[][3][3]){
    double C [2][3][3];
    double det[2];
    CofactorMatrix_3x3(A, C);
    Determinant_3x3(A, C, det);
    // Calculate B = adj(A)/det(A) = C^(T)/det(A)
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            B[0][i][j] = ((C[0][i][j]*det[0] - C[1][i][j]*det[1])/
                          (det[0]*det[0] + det[1]*det[1]));
            B[1][i][j] = -((C[0][i][j]*det[1] + C[1][i][j]*det[0])/
                           (det[0]*det[0] + det[1]*det[1]));
         }
    }
}

