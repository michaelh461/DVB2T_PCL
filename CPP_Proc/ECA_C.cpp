#include <math.h>
#include <algorithm>
#include "Cancellation.h"
#include "matrix_utils.h"

#define real_prod(a, b, c, d) (a*c - b*d)
#define comp_prod(a, b, c, d) (a*d + b*c)

double ** ECA_C(double ** RefSymbolData, double ** SurvSymbolData, int nCarriers, int nSymbols){
    double ** SurvSymbols_Cancelled = new double * [2];
    SurvSymbols_Cancelled[0] = new double[nSymbols*nCarriers];
    SurvSymbols_Cancelled[1] = new double[nSymbols*nCarriers];
    
    //Form doppler-shift matrix D
    double D[2][nSymbols];
    double Fd = 0.75;
    double T = 7/64e6;
    double Ts = nCarriers*T;
    const long double pi = 3.14159265358979323846264338328L;
   
    for(int k = 0; k < nCarriers; k++){
        // Reference channel carrier amplitudes (Q)
        double Q[2][nSymbols];
        std::copy(&RefSymbolData[0][0] + k*nSymbols, &RefSymbolData[0][0] + (k+1)*nSymbols, Q[0]);
        std::copy(&RefSymbolData[1][0] + k*nSymbols, &RefSymbolData[1][0] + (k+1)*nSymbols, Q[1]);
        // Surveillance channel carrier amplitudes (Y)
        double Y[2][nSymbols];
        std::copy(&SurvSymbolData[0][0] + k*nSymbols, &SurvSymbolData[0][0] + (k+1)*nSymbols, Y[0]);
        std::copy(&SurvSymbolData[1][0] + k*nSymbols, &SurvSymbolData[1][0] + (k+1)*nSymbols, Y[1]);
        
        double A = 0;
        double Ak = 0;
	for(int i = 0; i < nSymbols; i++){
		A = A + real_prod(Q[0][i], -Q[1][i], Q[0][i], Q[1][i]);
		Ak = Ak + comp_prod(Q[0][i], -Q[1][i], Q[0][i], Q[1][i]);
        }
        
        // Invert A matrix, C = Ak^(-1)
        A = 1/A;
        
        double B[2];
        B[0] = 0;
	B[1] = 1;
        //Form Bk matrix, Bk = X'*Y
        for(int i = 0; i < nSymbols; i++){
	   B[0] = B[0] + real_prod(Q[0][i], -Q[1][i], Y[0][i], Y[1][i]);
	   B[1] = B[1] + comp_prod(Q[0][i], -Q[1][i], Y[0][i], Y[1][i]);
        }
        
        //Calculate F = A*Bk = (C.')*Bk/
        double F[2];
        F[0] = A*B[0];
	F[1] = A*B[1];
        
        // Allocate memory for real and complex parts of cancelled data
        double Z[2][nSymbols];
       
        // Perform ECA on carrier
        for(int i = 0; i < nSymbols; i++){
                Z[0][i] = Y[0][i] - (Q[0][i]*F[0] - Q[1][i]*F[1]);
                Z[1][i] = Y[1][i] - (Q[0][i]*F[1] + Q[1][i]*F[0]);
        }
        
       // Copy data to be returned
        std::copy(Z[0], Z[0] + nSymbols, SurvSymbols_Cancelled[0] + k*nSymbols);
        std::copy(Z[1], Z[1] + nSymbols, SurvSymbols_Cancelled[1] + k*nSymbols);
    }
    return SurvSymbols_Cancelled;
}
