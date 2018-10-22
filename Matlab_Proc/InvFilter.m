function ARDMatrix = InvFilter(RefSymbolData,SurvSymbolData)

    % Transpose such that each row corresponds to a unique carrier
    RefSymbolData = RefSymbolData.';
    SurvSymbolData = SurvSymbolData.';
    
    % Calculate ARD by inverse filtering
    H = SurvSymbolData./RefSymbolData;  
    size(H)
  %  ARDMatrix = fftshift(fft2(H)).';
   ARDMatrix = fftshift(fft(fft(H).'));
end