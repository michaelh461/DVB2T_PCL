function [SurvSymbols_Cancelled] = ECA_C(RefSymbolData, SurvSymbolData, nSymbols, nCarriers)

% Z = SurvSymbolData_Cancelled
Z(nSymbols, nCarriers) = 0;

for k = 1:nCarriers
    % Get carrier amplitudes
    Q_k = RefSymbolData(:,k); 
    Y_k = SurvSymbolData(:,k);
    % Apply least squares regression on carrier
    Z_k = (eye(nSymbols) - Q_k*(Q_k' * Q_k)^(-1)*Q_k')*Y_k;
    Z(:, k) = Z_k;
end

% Return cancelled symbols
SurvSymbols_Cancelled = Z;


