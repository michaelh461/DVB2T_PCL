#ifndef CUDA_ECA_H
#define CUDA_ECA_H
#endif

__device__ void CofactorMatrix_3x3(double A[][3][3], double C[][3][3]);
__device__ void Determinant_3x3(double A [][3][3], double C[][3][3], double det[2]);
__device__ void InvertMatrix_3x3(double A [][3][3], double B[][3][3]);
