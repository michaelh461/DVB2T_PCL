#include "Cancellation.h"
#include "RCF.h"
#include <iostream>
#include <fstream>
#include "ARD_Makers.h"
#include <chrono>

int main(void){
    // File containing RCF data
	std::string fileName = "../frames/symbol_data.rcf";

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
    
    	const int nCarriers = 27841;
    	const int nSymbols = 60;
    
    	double ** RefSymbolData = new double * [2];
    	RefSymbolData[0] = new double[nSymbols*nCarriers];
    	RefSymbolData[1] = new double[nSymbols*nCarriers];
    	
    	double ** SurvSymbolData = new double * [2];
    	SurvSymbolData[0] = new double[nSymbols*nCarriers];
    	SurvSymbolData[1] = new double[nSymbols*nCarriers];
    
    	for(int i = 0; i < nSymbols*nCarriers; i++){
        	RefSymbolData[0][i] = RefData[2*i];
        	RefSymbolData[1][i] = RefData[2*i + 1];
        	SurvSymbolData[0][i] = SurvData[2*i];
        	SurvSymbolData[1][i] = SurvData[2*i + 1];
    	}
       
	auto startTime = std::chrono::high_resolution_clock::now(); 
    	double ** SurvSymbols_Cancelled = ECA_CD(RefSymbolData, SurvSymbolData, nCarriers, nSymbols);
	auto endTime = std::chrono::high_resolution_clock::now();
	auto totalTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
	std::cout << "ECA_C took: " << totalTime << " ms\n";	

	//Refotmat data for inverse filtering
	for(int i = 0; i < nCarriers; i++){
		for(int j = 0; j < nSymbols; j++){
			RefSymbolData[0][j*nCarriers + i] = RefData[2*(i*nSymbols + j)];
			RefSymbolData[1][j*nCarriers + i] = RefData[2*(i*nSymbols + j) + 1]; 
			SurvSymbolData[0][j*nCarriers + i] = /*SurvData[2*(i*nSymbols + j)];*/SurvSymbols_Cancelled[0][i*nSymbols + j];
			SurvSymbolData[1][j*nCarriers + i] = /*SurvData[2*(i*nSymbols + j)];*/SurvSymbols_Cancelled[1][i*nSymbols + j];
		}
	}
	
	startTime = std::chrono::high_resolution_clock::now();
	double ** ARDMatrix = inverse_filter(RefSymbolData, SurvSymbolData, nCarriers, nSymbols);
	endTime = std::chrono::high_resolution_clock::now();
	totalTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
	std::cout << "Inverse Filtering took: " << totalTime << " ms\n";

    	for(int i = 0; i < nSymbols*nCarriers; i++){
		RefData[2*i] = ARDMatrix[0][i];
		RefData[2*i + 1] = ARDMatrix[1][i]; 
        	SurvData[2*i] = SurvSymbols_Cancelled[0][i];
        	SurvData[2*i + 1] = SurvSymbols_Cancelled[1][i];
    	}

    
    	std::ofstream ofs;
   	 ofs.open("../frames/cancelled_rcf.rcf");
    	ofs << (*oRCF);
    
    
    	return 0;
    
}
