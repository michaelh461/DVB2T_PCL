#ifndef ARD_H
#define ARD_H
#endif

double ** inverse_filter(double ** RefSymbolData, double ** SurvSymbolData, int nCarriers, int nSymbols);
void FX_Batches(double ** RefData, double ** SurvData, int nSamples, int nBatches, float nSampBatches, int nRangeBins, int nDopplerBins, double *** ARDMatrix);
