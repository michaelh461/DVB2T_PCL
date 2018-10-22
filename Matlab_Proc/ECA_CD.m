function [SurvSymbols_Cancelled] = ECA_CD(RefSymbolData, SurvSymbolData, nSymbols, nCarriers)

% Z = SurvSymbolData_Cancelled
Z(nSymbols, nCarriers) = 0;

% Doppler shift matrix for clutter suppression
fd = 0.75;
T = 7/64e6;
Ts = nCarriers*T;
D = exp(1j*2*pi*fd*(0:nSymbols-1)*Ts)';  
D = diag(D);

for k = 1:nCarriers
    % Get carrier amplitudes
    Q_k = RefSymbolData(:,k); 
    X_k = [D'*Q_k, Q_k, D*Q_k];
    Y_k = SurvSymbolData(:,k);
    % Apply least squares regression
    Z_k = (eye(nSymbols) - X_k*(X_k' * X_k)^(-1)*X_k')*Y_k;
    Z(:, k) = Z_k;
end

% Return cancelled symbols
SurvSymbols_Cancelled = Z;