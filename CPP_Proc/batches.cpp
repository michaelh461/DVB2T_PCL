#include "ARD_Makers.h"
//Pre-processor directives for calculating real and complex products
// (a + ib)*(c + id)
#define real_prod(a, b, c, d) (a*c - b*d)
#define comp_prod(a, b, c, d) (a*d + b*c)

double windowBlackman(double n, double N);

void FX_Batches(double ** RefData, double ** SurvData, int nSamples, int nBatches, float nSampBatches, int nRangeBins, int nDopplerBins, double *** ARDMatrix){
    	//Batches parameters
    	int batch_stride = ceil(nSamples/nBatches);
    	int nBatchSamples = (int)ceil(batch_stride*nSampBatches);
    	
	// Allocate memory for FFT outputs
   	double RefFFT[2][nBatchSamples];
    	double SurvFFT[2][nBatchSamples];
 	
	// Build correlation matrix
    	double *** CorrelationMatrix = new double ** [2];
    	CorrelationMatrix[0] = new double * [nBatches];
    	CorrelationMatrix[1] = new double * [nBatches];
    
    	//FFTW type definitions
    	fftw_complex * in, * out;
    	in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*nBatchSamples);
    	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*nBatchSamples);
    	// Plan FFTs
	fftw_plan p;
    	fftw_plan q;
    	fftw_plan w;
    	p = fftw_plan_dft_1d(nBatchSamples, in, out, FFTW_FORWARD, FFTW_ESTIMATE); 
    	q = fftw_plan_dft_1d(nBatchSamples, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    	
        // Batch iteration
    	for(int batchNo = 0; batchNo < nBatches; batchNo++){
        	CorrelationMatrix[0][batchNo] = new double [nRangeBins];
        	CorrelationMatrix[1][batchNo] = new double [nRangeBins];
        	// Batch indices
        	int batchStartSample = batchNo*batch_stride;
        	int batchStopSample = batchStartSample + nBatchSamples;
        	if(batchStopSample > nSamples){
           	 	batchStopSample = nSamples;
            		//Zero-Padding as necessary
            		for(int i = batchStopSample-batchStartSample; i < nBatchSamples; i++){
                		in[i][0] = 0;
                		in[i][1] = 0;
            		}
        	}	 
        	// Take FFT of RefData
        	for(int i = batchStartSample; i < batchStopSample; i++){
            		in[i-batchStartSample][0] = RefData[0][i];
            		in[i-batchStartSample][1] = RefData[1][i];
        	}
        	fftw_execute(p);
        	for(int i = 0; i < nBatchSamples; i++){
            		RefFFT[0][i] = out[i][0];
            		RefFFT[1][i] = out[i][1];
        	}	 
        	// Take FFT of SurvData
        	for(int i = batchStartSample; i < batchStopSample; i++){
            		in[i-batchStartSample][0] = SurvData[0][i];
            		in[i-batchStartSample][1] = SurvData[1][i];
        	}
        	fftw_execute(p);
        	for(int i = 0; i < nBatchSamples; i++){
            		SurvFFT[0][i] = out[i][0];
            		SurvFFT[1][i] = out[i][1];
        	}
        	// Build correlation matrix, C(batchNo,:) = fft(SurvFFT .* conj(RefFFT))
        	for(int i = 0; i < nBatchSamples; i++){
            		in[i][0] = real_prod(SurvFFT[0][i], SurvFFT[1][i], RefFFT[0][i], -RefFFT[1][i]);
            		in[i][1] = comp_prod(SurvFFT[0][i], SurvFFT[1][i], RefFFT[0][i], -RefFFT[1][i]);
        	}
        	fftw_execute(q);
        	for(int i = 0; i < nRangeBins; i++){
            		CorrelationMatrix[0][batchNo][i] = out[i][0]/nBatchSamples;
            		CorrelationMatrix[1][batchNo][i] = out[i][1]/nBatchSamples;
        	}
    	}
    	// Build ARD Matrix
    	fftw_free(in);
    	fftw_free(out);
    	in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*nBatches);
    	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*nBatches);
    	w = fftw_plan_dft_1d(nBatches, in, out, FFTW_FORWARD, FFTW_ESTIMATE); 
    	for(int i = 0; i < nRangeBins; i++){
        	for(int j = 0; j < nBatches; j++){
            		// in*pow(-1, j) is equivalent to fftshift after taking fft
            		in[j][0] = CorrelationMatrix[0][j][i]*windowBlackman(j, nBatches)*pow(-1,j);
            		in[j][1] = CorrelationMatrix[1][j][i]*windowBlackman(j, nBatches)*pow(-1,j);
        	}
        	fftw_execute(w);
        	for(int j = 0; j < 2*nDopplerBins + 1; j++){
            		ARDMatrix[0][i][j] = out[(nBatches/2 - nDopplerBins) + j][0];
            		ARDMatrix[1][i][j] = out[(nBatches/2 - nDopplerBins) + j][1];
        	}
    	}	
    	//Free data
    	fftw_free(in);
    	fftw_free(out);
    	fftw_destroy_plan(p);
   	fftw_destroy_plan(q);
    	fftw_destroy_plan(w);
 //   for(int i = 0; i < nBatches; i++){
  //      delete [] CorrelationMatrix[0][i];
  //      delete [] CorrelationMatrix[1][i];
 //   }
 //   delete [] CorrelationMatrix[0];
 //   delete [] CorrelationMatrix[1];
}

double windowBlackman(double n, double N){
    const long double pi = 3.14159265358979323846264338328L;
    double factor = (2*pi*n)/(N - 1);
    return 0.42 - 0.5*cos(factor) + 0.08*cos(2*factor);
}