function out=HDM_OFT_CameraSpectralResponse(OFT_In_IDTTaskData)

OFT_Env=HDM_OFT_InitEnvironment();

HDM_OFT_Utils.OFT_DispTitle('camera spectral response estimation');

[l_Pixel2WavelengthLookUp_HigherOrderPolynomBased, l_LineR, l_LineG, l_LineB] = HDM_OFT_LineCalibrationPn...
    (OFT_In_IDTTaskData.SpectralResponse_In_LineCalibrationSpectrum, OFT_In_IDTTaskData.SpectralResponse_In_LineCalibrationImage,...
    OFT_In_IDTTaskData.PreLinearisation_Out_LinCurve, OFT_In_IDTTaskData.SpectralResponse_In_LightCalibrationImage);

OFT_CameraResponse = HDM_OFT_LightCalibrationPn...
    (l_Pixel2WavelengthLookUp_HigherOrderPolynomBased,...
    OFT_In_IDTTaskData.SpectralResponse_In_LightCalibrationSpectrum, OFT_In_IDTTaskData.SpectralResponse_In_LightCalibrationImage,...
    OFT_In_IDTTaskData.PreLinearisation_Out_LinCurve, OFT_In_IDTTaskData.Device_In_Sensor, OFT_In_IDTTaskData.Device_In_FocalLength);

%% error correction

l_IntenityAgainstWavelength = HDM_OFT_SpectrumExportImport.ImportSpectrum(OFT_In_IDTTaskData.SpectralResponse_In_LineCalibrationSpectrum);

OFT_Spectrum_MeanOfRowsRed_nmSamples=...
    interp1(l_Pixel2WavelengthLookUp_HigherOrderPolynomBased, l_LineR(1,:),...
    min(l_IntenityAgainstWavelength(1,:)):l_IntenityAgainstWavelength(1,2)-l_IntenityAgainstWavelength(1,1):max(l_IntenityAgainstWavelength(1,:)),...
    'linear');
OFT_Spectrum_MeanOfRowsGreen_nmSamples=...
    interp1(l_Pixel2WavelengthLookUp_HigherOrderPolynomBased, l_LineG(1,:),...
    min(l_IntenityAgainstWavelength(1,:)):l_IntenityAgainstWavelength(1,2)-l_IntenityAgainstWavelength(1,1):max(l_IntenityAgainstWavelength(1,:)),...
    'linear');
OFT_Spectrum_MeanOfRowsBlue_nmSamples=...
    interp1(l_Pixel2WavelengthLookUp_HigherOrderPolynomBased, l_LineB(1,:),...
    min(l_IntenityAgainstWavelength(1,:)):l_IntenityAgainstWavelength(1,2)-l_IntenityAgainstWavelength(1,1):max(l_IntenityAgainstWavelength(1,:)),...
    'linear');

%% blue

l_refSpectrumWithResponseScalingB = OFT_CameraResponse(4,:) .* l_IntenityAgainstWavelength(2,:)./max(l_IntenityAgainstWavelength(2,:));
l_refSpectrumWithResponseScalingB = l_refSpectrumWithResponseScalingB  ./ max(l_refSpectrumWithResponseScalingB(1,:));

figure('Name','Line Image Channnels with Response Scaling B')
plot(l_IntenityAgainstWavelength(1,:), l_refSpectrumWithResponseScalingB, 'b--', ...
    l_IntenityAgainstWavelength(1,:), OFT_Spectrum_MeanOfRowsBlue_nmSamples(1,:)./max(OFT_Spectrum_MeanOfRowsBlue_nmSamples(1,:)), 'b-')

%% blue corr lut

[OFT_SpectrumImage_peak_valueResponseScaledB, OFT_SpectrumImage_peak_locationResponseScaledB] = ...
    HDM_OFT_findpeaks(l_refSpectrumWithResponseScalingB,...
    5,10,100,...
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
    5,10,100,...
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

l_indexLUT_B = uint16(OFT_SpectrumImage_peak_locationResponseScaledB(1, 2:4));
l_waveLengthResponseScaledB = l_IntenityAgainstWavelength(1, l_indexLUT_B);

l_indexLUT_B = uint16(OFT_SpectrumImage_peak_locationCameraB(1, 1:3));
l_waveLengthCameraB = l_IntenityAgainstWavelength(1, l_indexLUT_B);

l_waveLength_B = 0.5 * (l_waveLengthResponseScaledB + l_waveLengthCameraB);
l_ratioCam2Scaled_B = OFT_SpectrumImage_peak_valueCameraB(1, 1:3) ./ OFT_SpectrumImage_peak_valueResponseScaledB(1,2:4);

l_corrLUT_B = [l_waveLength_B; l_ratioCam2Scaled_B];
[l_Y,l_I] = sort(l_corrLUT_B(1, :));
l_corrLUT_B = l_corrLUT_B(:, l_I);

%l_corrLUT_B = horzcat(l_corrLUT_B(1:2, 1:2), l_corrLUT_B(1:2, 4:4)); 

l_corrP1_B = polyfit(l_corrLUT_B(1, :), l_corrLUT_B(2, :), 1);

l_corrFullLUT_B = polyval(l_corrP1_B, OFT_CameraResponse(1, :));

l_corrFullLUTWithIndex_B = [OFT_CameraResponse(1, :); l_corrFullLUT_B];

%% green

l_refSpectrumWithResponseScalingG = OFT_CameraResponse(3,:) .* l_IntenityAgainstWavelength(2,:)./max(l_IntenityAgainstWavelength(2,:));
l_refSpectrumWithResponseScalingG = l_refSpectrumWithResponseScalingG  ./ max(l_refSpectrumWithResponseScalingG(1,:));

figure('Name','Line Image Channnels with Response Scaling G')
plot(l_IntenityAgainstWavelength(1,:), l_refSpectrumWithResponseScalingG, 'g--', ...
    l_IntenityAgainstWavelength(1,:), OFT_Spectrum_MeanOfRowsGreen_nmSamples(1,:)./max(OFT_Spectrum_MeanOfRowsGreen_nmSamples(1,:)), 'g-')

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
        AverageToClosePeaks(l_waveLengthResponseScaledG, OFT_SpectrumImage_peak_valueResponseScaledG, 6);

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
        AverageToClosePeaks(l_waveLengthCameraG, OFT_SpectrumImage_peak_valueCameraG, 6);

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

l_corrP1_G = polyfit(l_corrLUT_G(1, :), l_corrLUT_G(2, :), 1);

l_corrFullLUT_G = polyval(l_corrP1_G, OFT_CameraResponse(1, :));

l_corrFullLUTWithIndex_G = [OFT_CameraResponse(1, :); l_corrFullLUT_G];

%% red

l_refSpectrumWithResponseScalingR = OFT_CameraResponse(2,:) .* l_IntenityAgainstWavelength(2,:)./max(l_IntenityAgainstWavelength(2,:));
l_refSpectrumWithResponseScalingR = l_refSpectrumWithResponseScalingR  ./ max(l_refSpectrumWithResponseScalingR(1,:));

figure('Name','Line Image Channnels with Response Scaling R')
plot(l_IntenityAgainstWavelength(1,:), l_refSpectrumWithResponseScalingR, 'r--', ...
    l_IntenityAgainstWavelength(1,:), OFT_Spectrum_MeanOfRowsRed_nmSamples(1,:)./max(OFT_Spectrum_MeanOfRowsRed_nmSamples(1,:)), 'r-')

%% red corr lut

[OFT_SpectrumImage_peak_valueResponseScaledR, OFT_SpectrumImage_peak_locationResponseScaledR] = ...
    HDM_OFT_findpeaks(l_refSpectrumWithResponseScalingR,...
    5,10,100,...
    0.1,...
    false);

[l_peaklocs, l_sorter] = sort(OFT_SpectrumImage_peak_locationResponseScaledR);
l_amps = zeros(size(OFT_SpectrumImage_peak_valueResponseScaledR));
l_npeaks = length(l_peaklocs);
for cur = 1:l_npeaks,
    l_amps(:,cur) = OFT_SpectrumImage_peak_valueResponseScaledR(:, l_sorter(cur));
end;

OFT_SpectrumImage_peak_valueResponseScaledR = [l_amps(1, 1), l_amps(1, size(l_amps, 2))];
OFT_SpectrumImage_peak_locationResponseScaledR = [l_peaklocs(1, 1), l_peaklocs(1, size(l_peaklocs, 2))];

[OFT_SpectrumImage_peak_valueCameraR, OFT_SpectrumImage_peak_locationCameraR] = ...
    HDM_OFT_findpeaks(OFT_Spectrum_MeanOfRowsRed_nmSamples(1,:)./max(OFT_Spectrum_MeanOfRowsRed_nmSamples(1,:)),...
    5,10,100,...
    0.1,...
    false);

[l_peaklocs, l_sorter] = sort(OFT_SpectrumImage_peak_locationCameraR);
l_amps = zeros(size(OFT_SpectrumImage_peak_valueCameraR));
l_npeaks = length(l_peaklocs);
for cur = 1:l_npeaks,
    l_amps(:,cur) = OFT_SpectrumImage_peak_valueCameraR(:, l_sorter(cur));
end;

OFT_SpectrumImage_peak_valueCameraR = [l_amps(1, 1), l_amps(1, size(l_amps, 2))];
OFT_SpectrumImage_peak_locationCameraR = [l_peaklocs(1, 1), l_peaklocs(1, size(l_peaklocs, 2))];

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

l_corrP1_R = polyfit(l_corrLUT_R(1, :), l_corrLUT_R(2, :), 1);

l_corrFullLUT_R = polyval(l_corrP1_R, OFT_CameraResponse(1, :));

l_corrFullLUTWithIndex_R = [OFT_CameraResponse(1, :); l_corrFullLUT_R];

%% mean LUT

l_corrP1_Mean = (l_corrP1_R + l_corrP1_G + l_corrP1_B) / 3;
l_corrFullLUT_Mean = polyval(l_corrP1_Mean, OFT_CameraResponse(1, :));

l_corrFullLUTWithIndex_Mean = [OFT_CameraResponse(1, :); l_corrFullLUT_Mean];


%% apply corr LUT

out = OFT_CameraResponse;

out = [OFT_CameraResponse(1, :); ...
    l_corrFullLUTWithIndex_G(2, :) .* OFT_CameraResponse(2, :); ...
    l_corrFullLUTWithIndex_G(2, :) .* OFT_CameraResponse(3, :); ...
    l_corrFullLUTWithIndex_G(2, :) .* OFT_CameraResponse(4, :)];

out = [OFT_CameraResponse(1, :); ...
    out(2, :) / max(out(3, :)); ...
    out(3, :) / max(out(3, :)); ...
    out(4, :) / max(out(3, :))];

%% remove stray light error

l_index400nm = find (OFT_CameraResponse(1, :) == 400);

out(2, 1 : l_index400nm) = 0;
out(3, 1 : l_index400nm) = 0;
out(4, 1 : l_index400nm) = 0;

%% view

figure('Name','Corrected')
plot(OFT_CameraResponse(1, :), OFT_CameraResponse(2, :), 'r--', ...
    OFT_CameraResponse(1, :), OFT_CameraResponse(3, :), 'g--', ...
    OFT_CameraResponse(1, :), OFT_CameraResponse(4, :), 'b--', ...
    OFT_CameraResponse(1, :), out(2, :), 'r-', ...
    OFT_CameraResponse(1, :), out(3, :), 'g-', ...
    OFT_CameraResponse(1, :), out(4, :), 'b-');

HDM_OFT_Utils.OFT_DispTitle('camera spectral response estimation succesfully finished');

end

function [o_Peak_location, o_Peak_value] = AverageToClosePeaks(i_Peak_location, i_Peak_value, i_width)

l_PeakMatrix = [i_Peak_location; i_Peak_value];
[l_Y,l_I] = sort(l_PeakMatrix(1, :));
l_B = l_PeakMatrix(:, l_I); 

l_ReferencePeaksToCloseWaveLengthsFiltered = [];
l_ReferencePeaksToCloseWaveLengthsValuesFiltered = [];

l_peakIndex = 1;

while l_peakIndex <= (size(l_B, 2) - 1)
  
    if (abs(l_B(1, l_peakIndex + 1) - l_B(1, l_peakIndex)) > i_width)
    
        l_ReferencePeaksToCloseWaveLengthsFiltered = [l_ReferencePeaksToCloseWaveLengthsFiltered, l_B(1, l_peakIndex)];
        l_ReferencePeaksToCloseWaveLengthsValuesFiltered = [l_ReferencePeaksToCloseWaveLengthsValuesFiltered, l_B(2, l_peakIndex)];
        
        l_peakIndex = l_peakIndex + 1;
        
    elseif(abs(l_B(2, l_peakIndex) - l_B(2, l_peakIndex + 1)) > 0.5)
        
        l_ReferencePeaksToCloseWaveLengthsFiltered = [l_ReferencePeaksToCloseWaveLengthsFiltered, l_B(1, l_peakIndex)];
        l_ReferencePeaksToCloseWaveLengthsValuesFiltered = [l_ReferencePeaksToCloseWaveLengthsValuesFiltered, l_B(2, l_peakIndex)];
        
        l_peakIndex = l_peakIndex + 1;
        
    else    
        
        l_ReferencePeaksToCloseWaveLengthsFiltered = [l_ReferencePeaksToCloseWaveLengthsFiltered, (l_B(1, l_peakIndex) + l_B(1, l_peakIndex + 1)) / 2.0];
        l_ReferencePeaksToCloseWaveLengthsValuesFiltered = [l_ReferencePeaksToCloseWaveLengthsValuesFiltered, (l_B(2, l_peakIndex) + l_B(2, l_peakIndex + 1)) / 2.0];
        
        l_peakIndex = l_peakIndex + 2;
        
    end
    
end

if (abs(l_B(1, size(l_B, 2)) - l_B(1, size(l_B, 2) - 1)) > i_width)
    
    l_ReferencePeaksToCloseWaveLengthsFiltered = [l_ReferencePeaksToCloseWaveLengthsFiltered, l_B(1, size(l_B, 2))];
    l_ReferencePeaksToCloseWaveLengthsValuesFiltered = [l_ReferencePeaksToCloseWaveLengthsValuesFiltered, l_B(2, size(l_B, 2))];
    
elseif(abs(l_B(2, size(l_B, 2)) - l_B(2, size(l_B, 2) - 1)) > 0.5)
    
    l_ReferencePeaksToCloseWaveLengthsFiltered = [l_ReferencePeaksToCloseWaveLengthsFiltered, l_B(1, size(l_B, 2))];
    l_ReferencePeaksToCloseWaveLengthsValuesFiltered = [l_ReferencePeaksToCloseWaveLengthsValuesFiltered, l_B(2, size(l_B, 2))];

end

l_PeakMatrixToClosePeaksAveraged = [l_ReferencePeaksToCloseWaveLengthsValuesFiltered; l_ReferencePeaksToCloseWaveLengthsFiltered];
[l_Yout,l_Iout] = sort(l_PeakMatrixToClosePeaksAveraged(1, :), 'descend');
l_Bout = l_PeakMatrixToClosePeaksAveraged(:,l_Iout); 

o_Peak_value = l_Bout(1, :);
o_Peak_location = l_Bout(2, :);

[l_peaklocs, l_sorter] = sort(o_Peak_location);
l_amps = zeros(size(o_Peak_value));
l_npeaks = length(l_peaklocs);
for cur = 1:l_npeaks,
    l_amps(:,cur) = o_Peak_value(:, l_sorter(cur));
end;

o_Peak_value = l_amps;
o_Peak_location = l_peaklocs;

end
