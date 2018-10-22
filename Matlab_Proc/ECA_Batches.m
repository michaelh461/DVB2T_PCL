function [ SurvData_Cancelled ] = ECA_Batches(RefData_demod, SurvData_demod, nRangeBins, nDopplerBins, nBatches )

nSamples = max(size(RefData_demod));
nBatchSamples = floor((nSamples/nBatches));

X(nBatchSamples, nRangeBins) = 0;
SurvData_Cancelled(nSamples, 1) = 0;

for i = 1:nBatches
    batchStartSample = (i - 1)*nBatchSamples + 1;
    batchStopSample = i*nBatchSamples;
    if(batchStopSample > nSamples)
        batchStopSample = nSamples;
    end
    RefBatch = RefData_demod(batchStartSample : batchStopSample);
    SurvBatch = SurvData_demod(batchStartSample : batchStopSample);
    for k = 1:nRangeBins
        X(k:nBatchSamples,k)  = RefBatch(1:nBatchSamples - k + 1);
    end
    alpha = (X'*X)^(-1)*X'*SurvBatch;
    SurvData_Cancelled(batchStartSample : batchStopSample) = SurvBatch - X*alpha;
end



