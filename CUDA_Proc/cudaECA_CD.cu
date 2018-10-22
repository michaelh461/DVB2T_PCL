#include "RCF.h"
#include "cudaECA.cuh"
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include "ARD_Makers.h"

#define real_prod(a, b, c, d) (a*c - b*d)
#define comp_prod(a, b, c, d) (a*d + b*c)

using namespace std;

__global__ void ECA_CD(double * RefSymbolData, double * SurvSymbolData, double * D, double * SurvSymbols_Cancelled){
   
    const int nSymbols = 60;
    const int nCarriers = 27841;

    int index = blockIdx.x*blockDim.x + threadIdx.x;
    int stride = (int)nCarriers/(blockDim.x * gridDim.x);
 
    int start = index*stride;
    int stop = start + stride;
    if(stop > nCarriers){
        stop = nCarriers;
    }

    for(int k = start; k < stop; k++){
        // Reference channel carrier amplitudes (Q)
        double Q[2][nSymbols];
        // Surveillance channel carrier amplitudes (Y)
        double Y[2][nSymbols];
        for(int i = 0; i < nSymbols; i++){
            Q[0][i] = RefSymbolData[2*k*nSymbols + 2*i];
            Q[1][i] = RefSymbolData[2*k*nSymbols + 2*i + 1];
            Y[0][i] = SurvSymbolData[2*k*nSymbols + 2*i];
            Y[1][i] = SurvSymbolData[2*k*nSymbols + 2*i + 1];
        }
        
        // Form X matrix (clutter subspace matrix)
        double X[3][2][nSymbols];
        
        //X[0] = D' * Q
        for(int i = 0; i < nSymbols; i++){
            X[0][0][i] = real_prod(D[2*i], -D[2*i + 1], Q[0][i], Q[1][i]);
            X[0][1][i] = comp_prod(D[2*i], -D[2*i + 1], Q[0][i], Q[1][i]);
        }
        
        //X[1] = Q
        for(int i = 0; i < nSymbols; i++){
            X[1][0][i] =  Q[0][i];
            X[1][1][i] =  Q[1][i];
        }
    
        //X[2] = D * Q
        for(int i = 0; i < nSymbols; i++){
            X[2][0][i] = real_prod(D[2*i], D[2*i + 1], Q[0][i], Q[1][i]);
            X[2][1][i] = comp_prod(D[2*i], D[2*i + 1], Q[0][i], Q[1][i]);
        }    
        
        //Form A matrix; A = X'*X   
        double A[2][3][3];
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                A[0][i][j] = 0;
                A[1][i][j] = 0; 
                for(int n = 0; n < nSymbols; n++){
                    A[0][i][j] = A[0][i][j] + real_prod(X[j][0][n], -X[j][1][n], X[i][0][n], X[i][1][n]);
		    if(i == j){
			A[1][i][j] = 0;
		    }
		    else{
                        A[1][i][j] = A[1][i][j] + comp_prod(X[j][0][n], -X[j][1][n], X[i][0][n], X[i][1][n]);
		   }
                }                 
            }
        }
        
        // Invert A matrix, C = A^(-1)
        double C[2][3][3];
        InvertMatrix_3x3(A, C);
        
        //Form B matrix, B = X'*Y
        double B[2][3];
        for(int i = 0; i < 3; i++){
            B[0][i] = 0; 
            B[1][i] = 0;
            for(int n = 0; n < nSymbols; n++){
                B[0][i] = B[0][i] + X[i][0][n]*Y[0][n] + X[i][1][n]*Y[1][n];
                B[1][i] = B[1][i] + X[i][0][n]*Y[1][n] - X[i][1][n]*Y[0][n];
            }
        }
        
        //Calculate F = A*B = (C.')*B
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
        for(int i = 0; i < nSymbols; i++){
            SurvSymbols_Cancelled[2*k*nSymbols + 2*i] = Z[0][i];
            SurvSymbols_Cancelled[2*k*nSymbols + 2*i + 1] = Z[1][i];
        }
    }    

}

int main(void){
	// File containing RCF data
	string fileName = "../frames/symbol_data.rcf";

	// Create RCF object and read file header
	cRCF * oRCF = new cRCF();
	oRCF->readHeader(fileName, true);

	// Read data from file
	uint64_t nSamples = oRCF->getNSamples();
	oRCF->readData(fileName, 0, nSamples);

	// Get pointers to reference and surveillance samples
	float * RefData = oRCF->getReferenceArrayFloatPointer();
	float * SurvData = oRCF->getSurveillanceArrayFloatPointer();

	cout << "RCF read complete" << endl;
    
    const int nCarriers = 27841;
    const int nSymbols = 60;
    
    double * RefSymbolData;
    cudaMallocManaged(&RefSymbolData, 2*nSymbols*nCarriers*sizeof(double));
    double * SurvSymbolData;
    cudaMallocManaged(&SurvSymbolData, 2*nSymbols*nCarriers*sizeof(double));
    
    for(int i = 0; i < nSamples; i++){
        RefSymbolData[i*2] = (double)RefData[2*i];
        RefSymbolData[i*2+1] = (double)RefData[2*i+1];
        SurvSymbolData[i*2] = (double)SurvData[2*i];
        SurvSymbolData[i*2+1] = (double)SurvData[2*i+1];
    }

    //Form doppler-shift matrix D
    double * D;
    cudaMallocManaged(&D, 2*nSymbols*sizeof(double));
    // Note: these must be parametized
    double Fd = 0.75;
    double T = 7/64e6;
    double Ts = nCarriers*T;
    const long double pi = 3.14159265358979323846264338328L;
    for(int i = 0; i < nSymbols; i++){
        D[2*i] = cos(2*pi*Fd*i*Ts);
        D[2*i + 1] = sin(2*pi*Fd*i*Ts);
    }

    double * Surv_Cancelled;
    cudaMallocManaged(&Surv_Cancelled, 2*nSymbols*nCarriers*sizeof(double));
    
	int numSMs;
	cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, 0);
	std::cout << numSMs << std::endl;
    // Execute ECA_CD cuda core
    int n = 16*numSMs;
    ECA_CD<<<6, 256>>>(RefSymbolData, SurvSymbolData, D, Surv_Cancelled);
    cudaDeviceSynchronize();

    for(int i = 0; i < nSamples; i++){
        SurvData[2*i] = Surv_Cancelled[i*2];
        SurvData[2*i+1] = Surv_Cancelled[i*2 + 1];
    }

        std::ofstream ofs;
         ofs.open("../frames/cancelled_rcf.rcf");
        ofs << (*oRCF);


        return 0;



//    cudaFree(RefSymbolData);
//    cudaFree(SurvSymbolData);
//    cudaFree(D);
//    cudaFree(Surv_Cancelled);
    
//    ofstream ofs;
//    ofs.open("../frames/new_rcf.rcf");
//    ofs << (*oRCF);

    return 0;
}

// Utilities

/*  CofactorMatrix_3x3:
 *  Returns the matrix of cofactors of the 3x3 COMPLEX matrix A
 *  A[0] = 3x3 real data 
 *  A[1] = 3x3 imag data
 *  i.e. returns C such that adj(A) = C^(T)
 */
__device__ void CofactorMatrix_3x3(double A[][3][3], double C[][3][3]){
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
__device__ void Determinant_3x3(double A[][3][3], double C[][3][3], double det[2]){
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
__device__ void InvertMatrix_3x3(double A[][3][3], double B[][3][3]){
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
