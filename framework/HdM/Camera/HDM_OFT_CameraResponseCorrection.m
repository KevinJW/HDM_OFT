function o_correctedCameraResponse = HDM_OFT_CameraResponseCorrection...
    (i_firstEstimatedCameraResponse, ...
    i_lineCalibrationSpectrum, ...
    l_LineR, l_LineG, l_LineB, ...
    i_Pixel2WavelengthLookUp_HigherOrderPolynomBased)

HDM_OFT_Utils.OFT_DispTitle('begin camera response correction');

l_IntenityAgainstWavelength = HDM_OFT_SpectrumExportImport.ImportSpectrum(i_lineCalibrationSpectrum);

[l_LineSpectrumPeakValues, l_LineSpectrumPeakLocations] = ...
    HDM_OFT_findpeaks(l_IntenityAgainstWavelength(2, :),...
    5,10,100,...
    0.1,...
    false);

[l_peaklocs, l_sorter] = sort(l_LineSpectrumPeakLocations);
l_amps4Peaks = zeros(size(l_LineSpectrumPeakValues));
l_npeaks = length(l_peaklocs);
for cur = 1:l_npeaks,
    l_amps4Peaks(:,cur) = l_LineSpectrumPeakValues(:, l_sorter(cur));
end;

l_waveLengths4Peaks = l_IntenityAgainstWavelength(1, uint16(l_peaklocs(1, :)));

l_peakLocationsInGreen = [];
l_peakAmpsInGreen = [];

for cur = 1 : size(l_waveLengths4Peaks, 2)
    
    if (l_waveLengths4Peaks(1, cur) > 500 && l_waveLengths4Peaks(1, cur) < 600)
        
        l_peakLocationsInGreen = [l_peakLocationsInGreen, l_waveLengths4Peaks(1, cur)];
        l_peakAmpsInGreen = [l_peakAmpsInGreen, l_amps4Peaks(1, cur)];
        
    end
    
end

l_maxLocation = find(l_peakAmpsInGreen == max(l_peakAmpsInGreen(:)));
l_wavelength4Max = l_peakLocationsInGreen(l_maxLocation);

OFT_Spectrum_MeanOfRowsRed_nmSamples=...
    interp1(i_Pixel2WavelengthLookUp_HigherOrderPolynomBased, l_LineR(1,:),...
    min(l_IntenityAgainstWavelength(1,:)):l_IntenityAgainstWavelength(1,2)-l_IntenityAgainstWavelength(1,1):max(l_IntenityAgainstWavelength(1,:)),...
    'linear');
OFT_Spectrum_MeanOfRowsGreen_nmSamples=...
    interp1(i_Pixel2WavelengthLookUp_HigherOrderPolynomBased, l_LineG(1,:),...
    min(l_IntenityAgainstWavelength(1,:)):l_IntenityAgainstWavelength(1,2)-l_IntenityAgainstWavelength(1,1):max(l_IntenityAgainstWavelength(1,:)),...
    'linear');
OFT_Spectrum_MeanOfRowsBlue_nmSamples=...
    interp1(i_Pixel2WavelengthLookUp_HigherOrderPolynomBased, l_LineB(1,:),...
    min(l_IntenityAgainstWavelength(1,:)):l_IntenityAgainstWavelength(1,2)-l_IntenityAgainstWavelength(1,1):max(l_IntenityAgainstWavelength(1,:)),...
    'linear');

%% blue

l_refSpectrumWithResponseScalingB = i_firstEstimatedCameraResponse(4,:) .* l_IntenityAgainstWavelength(2,:)./max(l_IntenityAgainstWavelength(2,:));
l_refSpectrumWithResponseScalingB = l_refSpectrumWithResponseScalingB  ./ max(l_refSpectrumWithResponseScalingB(1,:));

HDM_OFT_UI_PlotAndSave...
    ([l_IntenityAgainstWavelength(1,:); ...
    l_refSpectrumWithResponseScalingB; OFT_Spectrum_MeanOfRowsBlue_nmSamples(1,:)./max(OFT_Spectrum_MeanOfRowsBlue_nmSamples(1,:))], ...
    'Line Image Blue Channnel', 'Wavelength (nm)' , 'Normalized Response', ...
    {'Spectrum Scaled by First Response' 'Spectrum from Image (Expectation)'})

%% blue corr lut

[OFT_SpectrumImage_peak_valueResponseScaledB, OFT_SpectrumImage_peak_locationResponseScaledB] = ...
    HDM_OFT_findpeaks(l_refSpectrumWithResponseScalingB,...
    3,10,100,...
    0.1,...
    false);

[l_peaklocs, l_sorter] = sort(OFT_SpectrumImage_peak_locationResponseScaledB);
l_amps = zeros(size(OFT_SpectrumImage_peak_valueResponseScaledB));
l_npeaks = length(l_peaklocs);
for cur = 1:l_npeaks,
    l_amps(:,cur) = OFT_SpectrumImage_peak_valueResponseScaledB(:, l_sorter(cur));
end;

OFT_SpectrumImage_peak_valueResponseScaledB = l_amps;
OFT_SpectrumImage_peak_locationResponseScaledB = l_peaklocs;

[OFT_SpectrumImage_peak_valueCameraB, OFT_SpectrumImage_peak_locationCameraB] = ...
    HDM_OFT_findpeaks(OFT_Spectrum_MeanOfRowsBlue_nmSamples(1,:)./max(OFT_Spectrum_MeanOfRowsBlue_nmSamples(1,:)),...
    3,10,100,...
    0.1,...
    false);

[l_peaklocs, l_sorter] = sort(OFT_SpectrumImage_peak_locationCameraB);
l_amps = zeros(size(OFT_SpectrumImage_peak_valueCameraB));
l_npeaks = length(l_peaklocs);
for cur = 1:l_npeaks,
    l_amps(:,cur) = OFT_SpectrumImage_peak_valueCameraB(:, l_sorter(cur));
end;

OFT_SpectrumImage_peak_valueCameraB = l_amps;
OFT_SpectrumImage_peak_locationCameraB = l_peaklocs;

l_indexLUT_B = uint16(OFT_SpectrumImage_peak_locationResponseScaledB(1, :));
l_waveLengthResponseScaledB = l_IntenityAgainstWavelength(1, l_indexLUT_B);

l_indexLUT_B = uint16(OFT_SpectrumImage_peak_locationCameraB(1, :));
l_waveLengthCameraB = l_IntenityAgainstWavelength(1, l_indexLUT_B);

l_waveLength_B = 0.5 * (l_waveLengthResponseScaledB + l_waveLengthCameraB);
l_ratioCam2Scaled_B = OFT_SpectrumImage_peak_valueCameraB(1, :) ./ OFT_SpectrumImage_peak_valueResponseScaledB(1,:);

l_corrLUT_B = [l_waveLength_B; l_ratioCam2Scaled_B];
[l_Y,l_I] = sort(l_corrLUT_B(1, :));
l_corrLUT_B = l_corrLUT_B(:, l_I);

%l_corrLUT_B = horzcat(l_corrLUT_B(1:2, 1:2), l_corrLUT_B(1:2, 4:4)); 

[l_c, l_index4ClosestGreenPeak] = min(abs(l_corrLUT_B(1, :) - l_wavelength4Max));

l_corrLUT_B(2, :) = l_corrLUT_B(2, :) / l_corrLUT_B(2, l_index4ClosestGreenPeak);

l_corrP1_B = polyfit(l_corrLUT_B(1, l_index4ClosestGreenPeak - 1 : l_index4ClosestGreenPeak), ...
    l_corrLUT_B(2, l_index4ClosestGreenPeak - 1 : l_index4ClosestGreenPeak), 1);

l_corrFullLUT_B = polyval(l_corrP1_B, i_firstEstimatedCameraResponse(1, :));

l_corrFullLUTWithIndex_B = [i_firstEstimatedCameraResponse(1, :); l_corrFullLUT_B];

%% green

l_refSpectrumWithResponseScalingG = i_firstEstimatedCameraResponse(3,:) .* l_IntenityAgainstWavelength(2,:)./max(l_IntenityAgainstWavelength(2,:));
l_refSpectrumWithResponseScalingG = l_refSpectrumWithResponseScalingG  ./ max(l_refSpectrumWithResponseScalingG(1,:));

HDM_OFT_UI_PlotAndSave...
    ([l_IntenityAgainstWavelength(1,:); ...
    l_refSpectrumWithResponseScalingG; OFT_Spectrum_MeanOfRowsGreen_nmSamples(1,:)./max(OFT_Spectrum_MeanOfRowsGreen_nmSamples(1,:))], ...
    'Line Image Green Channnel', 'Wavelength (nm)' , 'Normalized Response', ...
    {'Spectrum Scaled by First Response' 'Spectrum from Image (Expectation)'})

%% green corr lut

[OFT_SpectrumImage_peak_valueResponseScaledG, OFT_SpectrumImage_peak_locationResponseScaledG] = ...
    HDM_OFT_findpeaks(l_refSpectrumWithResponseScalingG,...
    5,10,100,...
    0.1,...
    false);

[l_peaklocs, l_sorter] = sort(OFT_SpectrumImage_peak_locationResponseScaledG);
l_amps = zeros(size(OFT_SpectrumImage_peak_valueResponseScaledG));
l_npeaks = length(l_peaklocs);
for cur = 1:l_npeaks,
    l_amps(:,cur) = OFT_SpectrumImage_peak_valueResponseScaledG(:, l_sorter(cur));
end;

OFT_SpectrumImage_peak_valueResponseScaledG = l_amps;
OFT_SpectrumImage_peak_locationResponseScaledG = l_peaklocs;

[OFT_SpectrumImage_peak_valueCameraG, OFT_SpectrumImage_peak_locationCameraG] = ...
    HDM_OFT_findpeaks(OFT_Spectrum_MeanOfRowsGreen_nmSamples(1,:)./max(OFT_Spectrum_MeanOfRowsGreen_nmSamples(1,:)),...
    5,10,100,...
    0.1,...
    false);

[l_peaklocs, l_sorter] = sort(OFT_SpectrumImage_peak_locationCameraG);
l_amps = zeros(size(OFT_SpectrumImage_peak_valueCameraG));
l_npeaks = length(l_peaklocs);
for cur = 1:l_npeaks,
    l_amps(:,cur) = OFT_SpectrumImage_peak_valueCameraG(:, l_sorter(cur));
end;

OFT_SpectrumImage_peak_valueCameraG = l_amps;
OFT_SpectrumImage_peak_locationCameraG = l_peaklocs;

l_indexLUT_G = uint16(OFT_SpectrumImage_peak_locationResponseScaledG(1,:));
l_waveLengthResponseScaledG = l_IntenityAgainstWavelength(1, l_indexLUT_G);

l_indexLUT_G = uint16(OFT_SpectrumImage_peak_locationCameraG(1,:));
l_waveLengthCameraG = l_IntenityAgainstWavelength(1, l_indexLUT_G);

% ignore peaks <450 + in 550/600

l_waveLengthResponseScaledGFiltered = [];
OFT_SpectrumImage_peak_valueResponseScaledGFiltered = [];
for cur = 1 : size(l_waveLengthResponseScaledG, 2)
    
    if (l_waveLengthResponseScaledG(cur) < 450 || (l_waveLengthResponseScaledG(cur) > 550 && l_waveLengthResponseScaledG(cur) < 600))
        continue;
    end
    
    l_waveLengthResponseScaledGFiltered = [l_waveLengthResponseScaledGFiltered, l_waveLengthResponseScaledG(cur)];
    OFT_SpectrumImage_peak_valueResponseScaledGFiltered = [OFT_SpectrumImage_peak_valueResponseScaledGFiltered, OFT_SpectrumImage_peak_valueResponseScaledG(cur)];
    
end

l_waveLengthResponseScaledG = l_waveLengthResponseScaledGFiltered;
OFT_SpectrumImage_peak_valueResponseScaledG = OFT_SpectrumImage_peak_valueResponseScaledGFiltered;

l_waveLengthCameraGFiltered = [];
OFT_SpectrumImage_peak_valueCameraGFiltered = [];
for cur = 1 : size(l_waveLengthCameraG, 2)
    
    if (l_waveLengthCameraG(cur) < 450 || (l_waveLengthCameraG(cur) > 550 && l_waveLengthCameraG(cur) < 600))
        continue;
    end
    
    l_waveLengthCameraGFiltered = [l_waveLengthCameraGFiltered, l_waveLengthCameraG(cur)];
    OFT_SpectrumImage_peak_valueCameraGFiltered = [OFT_SpectrumImage_peak_valueCameraGFiltered, OFT_SpectrumImage_peak_valueCameraG(cur)];
    
end

l_waveLengthCameraG = l_waveLengthCameraGFiltered;
OFT_SpectrumImage_peak_valueCameraG = OFT_SpectrumImage_peak_valueCameraGFiltered;

l_break = 0;

while l_break == 0
    
    [l_waveLengthResponseScaledG, OFT_SpectrumImage_peak_valueResponseScaledG] = ...
        HDM_OFT_AverageToClosePeaks(l_waveLengthResponseScaledG, OFT_SpectrumImage_peak_valueResponseScaledG, 6);

    l_ar = unique(uint16(l_waveLengthResponseScaledG));
       
    if size(l_ar, 2) < size(l_waveLengthResponseScaledG, 2)
        
        l_break = 0;
        
    else
        
        l_break = 1;
        
    end

end

l_break = 0;

while l_break == 0
    
    [l_waveLengthCameraG, OFT_SpectrumImage_peak_valueCameraG] = ...
        HDM_OFT_AverageToClosePeaks(l_waveLengthCameraG, OFT_SpectrumImage_peak_valueCameraG, 6);

    l_ar = unique(uint16(l_waveLengthCameraG));
       
    if size(l_ar, 2) < size(l_waveLengthCameraG, 2)
        
        l_break = 0;
        
    else
        
        l_break = 1;
        
    end

end

l_waveLength_G = 0.5 * (l_waveLengthResponseScaledG + l_waveLengthCameraG);
l_ratioCam2Scaled_G = OFT_SpectrumImage_peak_valueCameraG(1,:) ./ OFT_SpectrumImage_peak_valueResponseScaledG(1,:);

l_corrLUT_G = [l_waveLength_G; l_ratioCam2Scaled_G];
[l_Y,l_I] = sort(l_corrLUT_G(1, :));
l_corrLUT_G = l_corrLUT_G(:, l_I);

%l_corrLUT_G = horzcat(l_corrLUT_G(1:2, 1:2), l_corrLUT_G(1:2, 4:4)); 

[l_c, l_index4ClosestGreenPeak] = min(abs(l_corrLUT_G(1, :) - l_wavelength4Max));

l_corrLUT_G(2, :) = l_corrLUT_G(2, :) / l_corrLUT_G(2, l_index4ClosestGreenPeak);

l_corrP1_G = polyfit(l_corrLUT_G(1, :), l_corrLUT_G(2, :), 1);

l_corrFullLUT_G = polyval(l_corrP1_G, i_firstEstimatedCameraResponse(1, :));

l_corrFullLUTWithIndex_G = [i_firstEstimatedCameraResponse(1, :); l_corrFullLUT_G];

%% red

l_refSpectrumWithResponseScalingR = i_firstEstimatedCameraResponse(2,:) .* l_IntenityAgainstWavelength(2,:)./max(l_IntenityAgainstWavelength(2,:));
l_refSpectrumWithResponseScalingR = l_refSpectrumWithResponseScalingR  ./ max(l_refSpectrumWithResponseScalingR(1,:));

HDM_OFT_UI_PlotAndSave...
    ([l_IntenityAgainstWavelength(1,:); ...
    l_refSpectrumWithResponseScalingR; OFT_Spectrum_MeanOfRowsRed_nmSamples(1,:)./max(OFT_Spectrum_MeanOfRowsRed_nmSamples(1,:))], ...
    'Line Image Red Channnel', 'Wavelength (nm)' , 'Normalized Response', ...
    {'Spectrum Scaled by First Response' 'Spectrum from Image (Expectation)'})


%% red corr lut

[OFT_SpectrumImage_peak_valueResponseScaledR, OFT_SpectrumImage_peak_locationResponseScaledR] = ...
    HDM_OFT_findpeaks(l_refSpectrumWithResponseScalingR,...
    5,10,100,...
    0.1,...
    false);

l_maxPeak_valueResponseScaledR = OFT_SpectrumImage_peak_valueResponseScaledR(1, 1);
l_maxPeak_locationResponseScaledR = OFT_SpectrumImage_peak_locationResponseScaledR(1, 1);

[l_peaklocs, l_sorter] = sort(OFT_SpectrumImage_peak_locationResponseScaledR);
l_amps = zeros(size(OFT_SpectrumImage_peak_valueResponseScaledR));
l_npeaks = length(l_peaklocs);
for cur = 1:l_npeaks,
    l_amps(:,cur) = OFT_SpectrumImage_peak_valueResponseScaledR(:, l_sorter(cur));
end;

OFT_SpectrumImage_peak_valueResponseScaledR = [l_amps(1, 1), l_maxPeak_valueResponseScaledR];
OFT_SpectrumImage_peak_locationResponseScaledR = [l_peaklocs(1, 1), l_maxPeak_locationResponseScaledR];

[OFT_SpectrumImage_peak_valueCameraR, OFT_SpectrumImage_peak_locationCameraR] = ...
    HDM_OFT_findpeaks(OFT_Spectrum_MeanOfRowsRed_nmSamples(1,:)./max(OFT_Spectrum_MeanOfRowsRed_nmSamples(1,:)),...
    5,10,100,...
    0.1,...
    false);

l_maxPeak_valueCameraR = OFT_SpectrumImage_peak_valueCameraR(1, 1);
l_maxPeak_locationCameraR = OFT_SpectrumImage_peak_locationCameraR(1, 1);

[l_peaklocs, l_sorter] = sort(OFT_SpectrumImage_peak_locationCameraR);
l_amps = zeros(size(OFT_SpectrumImage_peak_valueCameraR));
l_npeaks = length(l_peaklocs);
for cur = 1:l_npeaks,
    l_amps(:,cur) = OFT_SpectrumImage_peak_valueCameraR(:, l_sorter(cur));
end;

OFT_SpectrumImage_peak_valueCameraR = [l_amps(1, 1), l_maxPeak_valueCameraR];
OFT_SpectrumImage_peak_locationCameraR = [l_peaklocs(1, 1), l_maxPeak_locationCameraR];

l_indexLUT_R = uint16(OFT_SpectrumImage_peak_locationResponseScaledR(1,:));
l_waveLengthResponseScaledR = l_IntenityAgainstWavelength(1, l_indexLUT_R);

l_indexLUT_R = uint16(OFT_SpectrumImage_peak_locationCameraR(1,:));
l_waveLengthCameraR = l_IntenityAgainstWavelength(1, l_indexLUT_R);

l_waveLength_R = 0.5 * (l_waveLengthResponseScaledR + l_waveLengthCameraR);
l_ratioCam2Scaled_R = OFT_SpectrumImage_peak_valueCameraR(1,:) ./ OFT_SpectrumImage_peak_valueResponseScaledR(1,:);

l_corrLUT_R = [l_waveLength_R; l_ratioCam2Scaled_R];
[l_Y,l_I] = sort(l_corrLUT_R(1, :));
l_corrLUT_R = l_corrLUT_R(:, l_I);

%l_corrLUT_R = horzcat(l_corrLUT_R(1:2, 1:2), l_corrLUT_R(1:2, 4:4)); 

[l_c, l_index4ClosestGreenPeak] = min(abs(l_corrLUT_R(1, :) - l_wavelength4Max));

l_corrLUT_R(2, :) = l_corrLUT_R(2, :) / l_corrLUT_R(2, l_index4ClosestGreenPeak);

l_corrP1_R = polyfit(l_corrLUT_R(1, :), l_corrLUT_R(2, :), 1);

l_corrFullLUT_R = polyval(l_corrP1_R, i_firstEstimatedCameraResponse(1, :));

l_corrFullLUTWithIndex_R = [i_firstEstimatedCameraResponse(1, :); l_corrFullLUT_R];

%% apply corr LUT

o_correctedCameraResponse = [i_firstEstimatedCameraResponse(1, :); ...
    l_corrFullLUTWithIndex_G(2, :) .* i_firstEstimatedCameraResponse(2, :); ...
    l_corrFullLUTWithIndex_G(2, :) .* i_firstEstimatedCameraResponse(3, :); ...
    l_corrFullLUTWithIndex_G(2, :) .* i_firstEstimatedCameraResponse(4, :)];

o_correctedCameraResponse = [i_firstEstimatedCameraResponse(1, :); ...
    o_correctedCameraResponse(2, :) / max(o_correctedCameraResponse(3, :)); ...
    o_correctedCameraResponse(3, :) / max(o_correctedCameraResponse(3, :)); ...
    o_correctedCameraResponse(4, :) / max(o_correctedCameraResponse(3, :))];

%% remove stray light error

l_index400nm = find(i_firstEstimatedCameraResponse(1, :) == 400);

o_correctedCameraResponse(2, 1 : l_index400nm) = 0;
o_correctedCameraResponse(3, 1 : l_index400nm) = 0;
o_correctedCameraResponse(4, 1 : l_index400nm) = 0;

%% view

l_plotArgs = {{i_firstEstimatedCameraResponse(1, :), i_firstEstimatedCameraResponse(2, :), 'r--'}; ...
    {i_firstEstimatedCameraResponse(1, :), i_firstEstimatedCameraResponse(3, :), 'g--'}; ...
    {i_firstEstimatedCameraResponse(1, :), i_firstEstimatedCameraResponse(4, :), 'b--'}; ...
    {i_firstEstimatedCameraResponse(1, :), o_correctedCameraResponse(2, :), 'r-'}; ...
    {i_firstEstimatedCameraResponse(1, :), o_correctedCameraResponse(3, :), 'g-'}; ...
    {i_firstEstimatedCameraResponse(1, :), o_correctedCameraResponse(4, :), 'b-'}};

HDM_OFT_UI_PlotAndSave...
    (l_plotArgs, ...
    'Corrected and Uncorrected Response', 'Wavelength (nm)' , 'Normalized Response', ...
    {'r uncorrected' 'g uncorrected' 'b uncorrected' 'r corrected' 'g corrected' 'b corrected'})


HDM_OFT_Utils.OFT_DispTitle('finish camera response correction');

end
