#include <math.h>
#include <algorithm>
#include "fftw3.h"
#include "ARD_Makers.h"

//Pre-processor directives for calculating real and complex products
// (a + ib)*(c + id)
#define real_prod(a, b, c, d) (a*c - b*d)
#define comp_prod(a, b, c, d) (a*d + b*c)

double ** inverse_filter(double ** RefSymbolData, double ** SurvSymbolData, int nCarriers, int nSymbols){
    fftw_complex * in, * out;
    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*nCarriers);
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*nCarriers);
    fftw_plan p;
    
    // Output data
    double ** ARDMatrix = new double * [2];
    ARDMatrix[0] = new double [nSymbols*nCarriers];
    ARDMatrix[1] = new double [nSymbols*nCarriers];
    
    double ** RefSymbol = new double * [2];
    RefSymbol[0] = new double [nCarriers];
    RefSymbol[1] = new double [nCarriers];
    
    double ** SurvSymbol = new double * [2];
    SurvSymbol[0] = new double [nCarriers];
    SurvSymbol[1] = new double [nCarriers];
    
    double *** RangeLines = new double ** [2];
    RangeLines[0] = new double * [nSymbols];
    RangeLines[1] = new double * [nSymbols];
    for(int i = 0; i < nSymbols; i++){
        RangeLines[0][i] = new double [nCarriers];
        RangeLines[1][i] = new double [nCarriers];
    }
    
    p = fftw_plan_dft_1d(nCarriers, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    
    for(int k = 0; k < nSymbols; k++){
        std::copy(&RefSymbolData[0][0] + k*nCarriers, &RefSymbolData[0][0] + (k+1)*nCarriers, RefSymbol[0]);
        std::copy(&RefSymbolData[1][0] + k*nCarriers, &RefSymbolData[1][0] + (k+1)*nCarriers, RefSymbol[1]);
        std::copy(&SurvSymbolData[0][0] + k*nCarriers, &SurvSymbolData[0][0] + (k+1)*nCarriers, SurvSymbol[0]);
        std::copy(&SurvSymbolData[1][0] + k*nCarriers, &SurvSymbolData[1][0] + (k+1)*nCarriers, SurvSymbol[1]);
        
        //Dot division
        for(int i = 0; i < nCarriers; i++){
            in[i][0] = (SurvSymbol[0][i]*RefSymbol[0][i] + SurvSymbol[1][i]*RefSymbol[1][i])/(RefSymbol[0][i]*RefSymbol[0][i] + RefSymbol[1][i]*RefSymbol[1][i]);
            in[i][1] = (SurvSymbol[1][i]*RefSymbol[0][i] - SurvSymbol[0][i]*RefSymbol[1][i])/(RefSymbol[0][i]*RefSymbol[0][i] + RefSymbol[1][i]*RefSymbol[1][i]);
        }
        
        fftw_execute(p);
        
        for(int i = 0; i < nCarriers; i++){
            RangeLines[0][k][i] = out[i][0];
            RangeLines[1][k][i] = out[i][1];
        }
    }
    
    fftw_free(in);
    fftw_free(out);
    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*nSymbols);
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*nSymbols);
    p = fftw_plan_dft_1d(nSymbols, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    
    for(int k = 0; k < nCarriers; k++){
        for(int i = 0; i < nSymbols; i++){
            in[i][0] = RangeLines[0][i][k];
            in[i][1] = RangeLines[1][i][k];
        }
        
        fftw_execute(p);
        
        for(int i = 0; i < nSymbols; i++){
            ARDMatrix[0][k*nSymbols + i] = out[i][0];
            ARDMatrix[1][k*nSymbols + i] = out[i][1];
        }
    }
    
    return ARDMatrix;
}
