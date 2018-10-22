#include <math.h>
#include <algorithm>
#include <fftw3.h>
#include "RCF.h"
#include <iostream>
#include <exception>
#include <stdio.h>
#include <chrono>
#include "ARD_Makers.h"

//Pre-processor directives for calculating real and complex products
// (a + ib)*(c + id)
#define real_prod(a, b, c, d) (a*c - b*d)
#define comp_prod(a, b, c, d) (a*d + b*c)

int main(void){
    	// File containing RCF data
	std::string fileName = "../frames/batches_data.rcf";

	// Create RCF object and read file header
	cRCF * oRCF = new cRCF();
	oRCF->readHeader(fileName, true);

	// Read data from file
	uint64_t nSamples = oRCF->getNSamples();
	oRCF->readData(fileName, 0, nSamples);

	// Get pointers to reference and surveillance samples
	float * RefData = oRCF->getReferenceArrayFloatPointer();
	float * SurvData = oRCF->getSurveillanceArrayFloatPointer();

	std::cout << "RCF read complete\n" << std::endl;
    
    	double ** RefData_demod = new double * [2];
    	RefData_demod[0] = new double[nSamples];
    	RefData_demod[1] = new double[nSamples];
    	
    	double ** SurvData_demod = new double * [2];
    	SurvData_demod[0] = new double[nSamples];
    	SurvData_demod[1] = new double[nSamples];
    
    	for(int i = 0; i < nSamples; i++){
        	RefData_demod[0][i] = RefData[2*i];
        	RefData_demod[1][i] = RefData[2*i + 1];
        	SurvData_demod[0][i] = SurvData[2*i];
        	SurvData_demod[1][i] = SurvData[2*i + 1];
   	}
    
    	//Batch Parameters
    	int nBatches = 256;
    	float nSampBatches = 1.1;
	int ARDMaxRange_m = 300000;    
	int ARDMaxDoppler_Hz = 100;

	// Constant parameters
	int TxToRefRxDistance_m = 12600;
	int c = 299792458;
	double Fs = (64e6)/7;

    	int batch_stride = ceil(nSamples/nBatches);
    	int nBatchSamples = ceil(batch_stride*nSampBatches);
    	int nRangeBins = nBatchSamples;//ceil((ARDMaxRange_m - TxToRefRxDistance_m)*Fs/c);
    	// Calculate Doppler Bins
    	int nDopplerBins = 30;// floor(ARDMaxDoppler_Hz/((Fs/batch_stride)/nBatches));
        
      	std::cout << "Executing Batches" << std::endl;
	std::cout << "nRangeBins: " << nRangeBins << std::endl;
	std::cout << "nDopplerBins: " << nDopplerBins << std::endl;

      	double *** ARDMatrix = new double ** [2];
      	ARDMatrix[0] = new double * [nRangeBins];
      	ARDMatrix[1] = new double * [nRangeBins];
      	for(int i = 0; i < nRangeBins; i++){
         	ARDMatrix[0][i] = new double[2*nDopplerBins + 1];
         	ARDMatrix[1][i] = new double[2*nDopplerBins + 1];
      	}

      	auto t1 = std::chrono::high_resolution_clock::now();
      	FX_Batches(RefData_demod, SurvData_demod, nSamples, nBatches, nSampBatches, nRangeBins, nDopplerBins, ARDMatrix);
      	auto t2 = std::chrono::high_resolution_clock::now();
     	std::cout << "FX Batches took: "
               	  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
                  << " ms\n";
    
    	for(int i = 0; i < nRangeBins; i++){
		for(int j = 0; j < 2*nDopplerBins + 1; j++){
        		SurvData[2*(i*(2*nDopplerBins + 1) + j)] = ARDMatrix[0][i][j];
        		SurvData[2*(i*(2*nDopplerBins + 1) + j) + 1] = ARDMatrix[1][i][j];
		}
    	}	
    
    	std::ofstream ofs;
   	ofs.open("../frames/new_rcf.rcf");
    	ofs << (*oRCF);
    
    	return 0;
}

