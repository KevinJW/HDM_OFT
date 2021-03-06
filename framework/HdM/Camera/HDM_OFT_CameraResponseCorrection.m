function o_correctedCameraResponse = HDM_OFT_CameraResponseCorrection...
    (i_firstEstimatedCameraResponse, ...
    i_lineCalibrationSpectrum, ...
    l_LineR, l_LineG, l_LineB, ...
    i_Pixel2WavelengthLookUp_HigherOrderPolynomBased, ...
    i_lightCalibrationSpectrum, ...
    i_LightR, i_LightG, i_LightB)

HDM_OFT_Utils.OFT_DispTitle('begin camera response correction');

l_IntenityAgainstWavelengthLight = HDM_OFT_SpectrumExportImport.ImportSpectrum(i_lightCalibrationSpectrum);

l_refSpectrumWithResponseScalingB = i_firstEstimatedCameraResponse(4,:) .* ...
                                    l_IntenityAgainstWavelengthLight(2,:);

l_refSpectrumWithResponseScalingG = i_firstEstimatedCameraResponse(3,:) .* ...
                                    l_IntenityAgainstWavelengthLight(2,:);

l_refSpectrumWithResponseScalingR = i_firstEstimatedCameraResponse(2,:) .* ...
                                    l_IntenityAgainstWavelengthLight(2,:);

l_Spectrum_MeanOfRowsBlue_nmSamples=...
    interp1(i_Pixel2WavelengthLookUp_HigherOrderPolynomBased, i_LightB(1,:),...
    min(l_IntenityAgainstWavelengthLight(1,:)):l_IntenityAgainstWavelengthLight(1,2) - ...
    l_IntenityAgainstWavelengthLight(1,1):max(l_IntenityAgainstWavelengthLight(1,:)),...
    'linear');

l_Spectrum_MeanOfRowsGreen_nmSamples=...
    interp1(i_Pixel2WavelengthLookUp_HigherOrderPolynomBased, i_LightG(1,:),...
    min(l_IntenityAgainstWavelengthLight(1,:)):l_IntenityAgainstWavelengthLight(1,2) - ...
    l_IntenityAgainstWavelengthLight(1,1):max(l_IntenityAgainstWavelengthLight(1,:)),...
    'linear');

l_Spectrum_MeanOfRowsRed_nmSamples=...
    interp1(i_Pixel2WavelengthLookUp_HigherOrderPolynomBased, i_LightR(1,:),...
    min(l_IntenityAgainstWavelengthLight(1,:)):l_IntenityAgainstWavelengthLight(1,2) - ...
    l_IntenityAgainstWavelengthLight(1,1):max(l_IntenityAgainstWavelengthLight(1,:)),...
    'linear');

%% view

HDM_OFT_UI_PlotAndSave...
    ([l_IntenityAgainstWavelengthLight(1,:); ...
    l_refSpectrumWithResponseScalingB ./ max(l_refSpectrumWithResponseScalingG); ...
    l_Spectrum_MeanOfRowsBlue_nmSamples(1,:) ./ max(l_Spectrum_MeanOfRowsGreen_nmSamples(1,:));
    l_refSpectrumWithResponseScalingG ./ max(l_refSpectrumWithResponseScalingG); ...
    l_Spectrum_MeanOfRowsGreen_nmSamples(1,:) ./ max(l_Spectrum_MeanOfRowsGreen_nmSamples(1,:));
    l_refSpectrumWithResponseScalingR ./ max(l_refSpectrumWithResponseScalingG); ...
    l_Spectrum_MeanOfRowsRed_nmSamples(1,:) ./ max(l_Spectrum_MeanOfRowsGreen_nmSamples(1,:))], ...
    'Light Image Channnels', 'Wavelength (nm)' , 'Normalized Response', ...
    {'Spectrum Scaled by First Response B' 'Spectrum from Image (Expectation) B' ...
    'Spectrum Scaled by First Response G' 'Spectrum from Image (Expectation) G'...
    'Spectrum Scaled by First Response R' 'Spectrum from Image (Expectation) R'})

%% get shift via crosscorrelation

[l_acor,l_lag] = xcorr(l_Spectrum_MeanOfRowsBlue_nmSamples(1,:) ./ max(l_Spectrum_MeanOfRowsGreen_nmSamples(1,:)), ...
                        l_refSpectrumWithResponseScalingB ./ max(l_refSpectrumWithResponseScalingG));
[~,l_I] = max(abs(l_acor));
l_lagDiffB = l_lag(l_I);

[l_acor,l_lag] = xcorr(l_Spectrum_MeanOfRowsGreen_nmSamples(1,:) ./ max(l_Spectrum_MeanOfRowsGreen_nmSamples(1,:)), ...
                        l_refSpectrumWithResponseScalingG ./ max(l_refSpectrumWithResponseScalingG));
[~,l_I] = max(abs(l_acor));
l_lagDiffG = l_lag(l_I);

[l_acor,l_lag] = xcorr(l_Spectrum_MeanOfRowsRed_nmSamples(1,:) ./ max(l_Spectrum_MeanOfRowsGreen_nmSamples(1,:)), ...
                        l_refSpectrumWithResponseScalingR ./ max(l_refSpectrumWithResponseScalingG));
[~,l_I] = max(abs(l_acor));
l_lagDiffR = l_lag(l_I);

%shift currently not used
l_lagDiffR = 0;
l_lagDiffG = 0;
l_lagDiffB = 0;

l_firstEstimatedCameraResponseShiftedR = circshift(i_firstEstimatedCameraResponse(2, :), l_lagDiffR, 2);
l_firstEstimatedCameraResponseShiftedG = circshift(i_firstEstimatedCameraResponse(3, :), l_lagDiffG, 2);
l_firstEstimatedCameraResponseShiftedB = circshift(i_firstEstimatedCameraResponse(4, :), l_lagDiffB, 2);

l_firstEstimatedCameraResponseShifted = ...
    [i_firstEstimatedCameraResponse(1, :); ...
    l_firstEstimatedCameraResponseShiftedR; l_firstEstimatedCameraResponseShiftedG; l_firstEstimatedCameraResponseShiftedB];

%% calc shift corrected

l_refSpectrumWithResponseScalingB = l_firstEstimatedCameraResponseShifted(4,:) .* ...
                                    l_IntenityAgainstWavelengthLight(2,:);

l_refSpectrumWithResponseScalingG = l_firstEstimatedCameraResponseShifted(3,:) .* ...
                                    l_IntenityAgainstWavelengthLight(2,:);

l_refSpectrumWithResponseScalingR = l_firstEstimatedCameraResponseShifted(2,:) .* ...
                                    l_IntenityAgainstWavelengthLight(2,:);

%% view

HDM_OFT_UI_PlotAndSave...
    ([l_IntenityAgainstWavelengthLight(1,:); ...
    l_refSpectrumWithResponseScalingB ./ max(l_refSpectrumWithResponseScalingG); ...
    l_Spectrum_MeanOfRowsBlue_nmSamples(1,:) ./ max(l_Spectrum_MeanOfRowsGreen_nmSamples(1,:));
    l_refSpectrumWithResponseScalingG ./ max(l_refSpectrumWithResponseScalingG); ...
    l_Spectrum_MeanOfRowsGreen_nmSamples(1,:) ./ max(l_Spectrum_MeanOfRowsGreen_nmSamples(1,:));
    l_refSpectrumWithResponseScalingR ./ max(l_refSpectrumWithResponseScalingG); ...
    l_Spectrum_MeanOfRowsRed_nmSamples(1,:) ./ max(l_Spectrum_MeanOfRowsGreen_nmSamples(1,:))], ...
    'Light Image Channnels Corrected by Shift', 'Wavelength (nm)' , 'Normalized Response', ...
    {'Spectrum Scaled by First Response B' 'Spectrum from Image (Expectation) B' ...
    'Spectrum Scaled by First Response G' 'Spectrum from Image (Expectation) G'...
    'Spectrum Scaled by First Response R' 'Spectrum from Image (Expectation) R'})


%% calculate lin correction by blue and red peak ratio for green normalized samples

l_curB = l_refSpectrumWithResponseScalingB;
l_refB = l_Spectrum_MeanOfRowsBlue_nmSamples(1,:);
l_blueRatio = max(max(l_refB)/max(l_curB));

l_curG = l_refSpectrumWithResponseScalingG;
l_refG = l_Spectrum_MeanOfRowsGreen_nmSamples(1,:);
l_greenRatio = max(max(l_refG)/max(l_curG));

l_curR = l_refSpectrumWithResponseScalingR;
l_refR = l_Spectrum_MeanOfRowsRed_nmSamples(1,:);
l_redRatio = max(max(l_refR)/max(l_curR));

[l_num l_idx] = max(l_refB(:));
[l_x l_y] = ind2sub(size(l_refB),l_idx);

l_lambdaB = i_firstEstimatedCameraResponse(1, l_y);

[l_num l_idx] = max(l_refG(:));
[l_x l_y] = ind2sub(size(l_refG),l_idx);

l_lambdaG = i_firstEstimatedCameraResponse(1, l_y);

[l_num l_idx] = max(l_refR(:));
[l_x l_y] = ind2sub(size(l_refR),l_idx);

l_lambdaR = i_firstEstimatedCameraResponse(1, l_y);

% l_blueRatio = 0.85 * l_blueRatio;
% l_redRatio = 1.09 * l_redRatio;

l_corrLUT = interp1([l_lambdaB; l_lambdaG; l_lambdaR], [l_blueRatio; l_greenRatio; l_redRatio], i_firstEstimatedCameraResponse(1, :)', ...
                    'linear', 'extrap');
l_corrLUT = [i_firstEstimatedCameraResponse(1, :); l_corrLUT'];

%%%%%

o_correctedCameraResponse = l_firstEstimatedCameraResponseShifted;

o_correctedCameraResponse(4,:) = o_correctedCameraResponse(4,:) .* l_corrLUT(2,:);
o_correctedCameraResponse(3,:) = o_correctedCameraResponse(3,:) .* l_corrLUT(2,:);
o_correctedCameraResponse(2,:) = o_correctedCameraResponse(2,:) .* l_corrLUT(2,:);

l_norm = max(o_correctedCameraResponse(3,:));

o_correctedCameraResponse(3,:) = o_correctedCameraResponse(3,:) ./ l_norm;

o_correctedCameraResponse(2,:) = o_correctedCameraResponse(2,:) ./ l_norm;
o_correctedCameraResponse(4,:) = o_correctedCameraResponse(4,:) ./ l_norm;

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
    {'r uncorrected' 'g uncorrected' 'b uncorrected' 'r corrected' 'g corrected' 'b corrected'});

%% check again

l_refSpectrumWithResponseScalingB = o_correctedCameraResponse(4,:) .* ...
                                    l_IntenityAgainstWavelengthLight(2,:);

l_refSpectrumWithResponseScalingG = o_correctedCameraResponse(3,:) .* ...
                                    l_IntenityAgainstWavelengthLight(2,:);

l_refSpectrumWithResponseScalingR = o_correctedCameraResponse(2,:) .* ...
                                    l_IntenityAgainstWavelengthLight(2,:);                                

%% view

HDM_OFT_UI_PlotAndSave...
    ([l_IntenityAgainstWavelengthLight(1,:); ...
    l_refSpectrumWithResponseScalingB ./ max(l_refSpectrumWithResponseScalingG); ...
    l_Spectrum_MeanOfRowsBlue_nmSamples(1,:) ./ max(l_Spectrum_MeanOfRowsGreen_nmSamples(1,:));
    l_refSpectrumWithResponseScalingG ./ max(l_refSpectrumWithResponseScalingG); ...
    l_Spectrum_MeanOfRowsGreen_nmSamples(1,:) ./ max(l_Spectrum_MeanOfRowsGreen_nmSamples(1,:));
    l_refSpectrumWithResponseScalingR ./ max(l_refSpectrumWithResponseScalingG); ...
    l_Spectrum_MeanOfRowsRed_nmSamples(1,:) ./ max(l_Spectrum_MeanOfRowsGreen_nmSamples(1,:))], ...
    'Light Image Channnels By Corrected Response', 'Wavelength (nm)' , 'Normalized Response', ...
    {'Spectrum Scaled by Corrected Response B' 'Spectrum from Image (Expectation) B' ...
    'Spectrum Scaled by Corrected Response G' 'Spectrum from Image (Expectation) G'...
    'Spectrum Scaled by Corrected Response R' 'Spectrum from Image (Expectation) R'})

HDM_OFT_Utils.OFT_DispTitle('finish camera response correction');

end
