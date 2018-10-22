#include <math.h>
#include <algorithm>
#include "Cancellation.h"
#include "matrix_utils.h"

#define real_prod(a, b, c, d) (a*c - b*d)
#define comp_prod(a, b, c, d) (a*d + b*c)

double ** ECA_CD(double ** RefSymbolData, double ** SurvSymbolData, int nCarriers, int nSymbols){
    double ** SurvSymbols_Cancelled = new double * [2];
    SurvSymbols_Cancelled[0] = new double[nSymbols*nCarriers];
    SurvSymbols_Cancelled[1] = new double[nSymbols*nCarriers];
    
    //Form doppler-shift matrix D
    double D[2][nSymbols];
    double Fd = 0.75;
    double T = 7/64e6;
    double Ts = nCarriers*T;
    const long double pi = 3.14159265358979323846264338328L;
    for(int i = 0; i < nSymbols; i++){
        D[0][i] = cos(2*pi*Fd*i*Ts);
        D[1][i] = sin(2*pi*Fd*i*Ts);
    }
   
    for(int k = 0; k < nCarriers; k++){
        // Reference channel carrier amplitudes (Q)
        double Q[2][nSymbols];
        std::copy(&RefSymbolData[0][0] + k*nSymbols, &RefSymbolData[0][0] + (k+1)*nSymbols, Q[0]);
        std::copy(&RefSymbolData[1][0] + k*nSymbols, &RefSymbolData[1][0] + (k+1)*nSymbols, Q[1]);
        // Surveillance channel carrier amplitudes (Y)
        double Y[2][nSymbols];
        std::copy(&SurvSymbolData[0][0] + k*nSymbols, &SurvSymbolData[0][0] + (k+1)*nSymbols, Y[0]);
        std::copy(&SurvSymbolData[1][0] + k*nSymbols, &SurvSymbolData[1][0] + (k+1)*nSymbols, Y[1]);
        
        // Form X matrix (clutter subspace matrix)
        double X[3][2][nSymbols];
        
        //X[0] = D' * Q
        for(int i = 0; i < nSymbols; i++){
            X[0][0][i] = real_prod(D[0][i], -D[1][i], Q[0][i], Q[1][i]);
            X[0][1][i] = comp_prod(D[0][i], -D[1][i], Q[0][i], Q[1][i]);
        }
        
        //X[1] = Q
        std::copy(&RefSymbolData[0][0] + k*nSymbols, &RefSymbolData[0][0] + (k+1)*nSymbols, X[1][0]);
        std::copy(&RefSymbolData[1][0] + k*nSymbols, &RefSymbolData[1][0] + (k+1)*nSymbols, X[1][1]);
    
        //X[2] = D * Q
        for(int i = 0; i < nSymbols; i++){
            X[2][0][i] = real_prod(D[0][i], D[1][i], Q[0][i], Q[1][i]);
            X[2][1][i] = comp_prod(D[0][i], D[1][i], Q[0][i], Q[1][i]);
        }    
        
        //Form A matrix; A = X'*X   
        double A[2][3][3];
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                A[0][i][j] = 0;
                A[1][i][j] = 0; 
                for(int n = 0; n < nSymbols; n++){
                    A[0][i][j] = A[0][i][j] + real_prod(X[j][0][n], -X[j][1][n], X[i][0][n], X[i][1][n]);
		    // Eliminate redundant calculation
		    if(i == j){
                         A[1][i][j] = 0;
		    }
		    else{
		         A[1][i][j] = A[1][i][j] + comp_prod(X[j][0][n], -X[j][1][n], X[i][0][n], X[i][1][n]);
		    }                
		}                 
            }
        }
        
        // Invert A matrix, C = Ak^(-1)
        double C[2][3][3];
        InvertMatrix_3x3(A, C);
        
        //Form Bk matrix, Bk = X'*Y
        double B[2][3];
        for(int i = 0; i < 3; i++){
            B[0][i] = 0; 
            B[1][i] = 0;
            for(int n = 0; n < nSymbols; n++){
                B[0][i] = B[0][i] + X[i][0][n]*Y[0][n] + X[i][1][n]*Y[1][n];
                B[1][i] = B[1][i] + X[i][0][n]*Y[1][n] - X[i][1][n]*Y[0][n];
            }
        }
        
        //Calculate F = A*Bk = (C.')*Bk/
        double F[2][3];
        for(int i = 0; i < 3; i++){
            F[0][i] = 0; 
            F[1][i] = 0;
            for(int n = 0; n < 3; n++){
                F[0][i] = F[0][i] + C[0][n][i]*B[0][n] - C[1][n][i]*B[1][n];
                F[1][i] = F[1][i] + C[0][n][i]*B[1][n] + C[1][n][i]*B[0][n];
            }
        }
        
        // Allocate memory for real and complex parts of cancelled data
        double Z[2][nSymbols];
       
        // Perform ECA on carrier
        for(int i = 0; i < nSymbols; i++){
            Z[0][i] = Y[0][i];
            Z[1][i] = Y[1][i];
            for(int n = 0; n < 3; n++){
                Z[0][i] = Z[0][i] - (X[n][0][i]*F[0][n] - X[n][1][i]*F[1][n]);
                Z[1][i] = Z[1][i] - (X[n][0][i]*F[1][n] + X[n][1][i]*F[0][n]);
            }
        }
        
       // Copy data to be returned
        std::copy(Z[0], Z[0] + nSymbols, SurvSymbols_Cancelled[0] + k*nSymbols);
        std::copy(Z[1], Z[1] + nSymbols, SurvSymbols_Cancelled[1] + k*nSymbols);
    }
    return SurvSymbols_Cancelled;
}
