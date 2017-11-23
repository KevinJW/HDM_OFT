function [o_IDT_B, o_IDT_b, o_referenceWhite, o_resnormBEstimation, o_meanDeltaE2000, o_stdDevDeltaE2000]=HDM_OFT_IDT_MinimumError...
    (i_ObjectReflectances, ...
    i_CurrentObserver, ...
    i_PreLinearisationCurve4CurrentObserverImage, ...
    i_NeutralsCompensation, ...
    i_ErrorMinimizationDomain, ...
    i_ReferenceIlluminantSpectrum, ...
    i_CurrentIlluminantSpectrum, ...
    i_ReferenceDomain, ...
    i_ReferenceObserver, ...
    i_IndiceesOfObjectReflectances4CorrectMappedColour, ...
    i_indexOfReferenceWhite, ...
    i_ObjectReflectancesWeightings)

HDM_OFT_Utils.OFT_DispTitle('start matrix creation'); 

%%reference camera
HDM_OFT_Utils.OFT_DispSubTitle('setup reference camera');

%XYZ to RICD 
[l_M, l_w] = HDM_OFT_IDT_ReferenceCamera.GetDefinition(i_ReferenceDomain);

%% illuminat spectrum aquisition
HDM_OFT_Utils.OFT_DispSubTitle('illuminat spectrum aquisition');
OFT_ReferenceIlluminant_Spectrum_1nm_CIE31Range = HDM_OFT_GetIlluminantSpectrum(i_ReferenceIlluminantSpectrum);
OFT_CurrentIlluminant_Spectrum_1nm_CIE31Range = HDM_OFT_GetIlluminantSpectrum(i_CurrentIlluminantSpectrum);

%%patches spectrum aquisition
%10 nm resolution
HDM_OFT_Utils.OFT_DispSubTitle('read patch spectra');
    
OFT_PatchSet_SpectralCurve = HDM_OFT_PatchSet.GetPatchSpectra(i_ObjectReflectances);

par4gOFT_w = l_w;

%% visualize spectral reflectance subset for colorchecker

HDM_OFT_UI_PlotAndSave([OFT_PatchSet_SpectralCurve(1,:); OFT_PatchSet_SpectralCurve(2,:); OFT_PatchSet_SpectralCurve(7,:); OFT_PatchSet_SpectralCurve(19,:)], ...
        'Patch Reflectances for Selected Improved Colours', 'Wavelength (nm)', 'Normalized Reflectance', ...
        {'Dark Skin','Orange','Cyan'});

%% visualize SPD of used illuminant, CCT related spectrum and target illuminant

l_CIEStandardObserver_SpectralCurves = HDM_OFT_CIEStandard.GetStandardObserverCurves(i_ReferenceObserver);

l_XusedIllum = trapz(l_CIEStandardObserver_SpectralCurves(2,:) .* OFT_CurrentIlluminant_Spectrum_1nm_CIE31Range(2,:));
l_YusedIllum = trapz(l_CIEStandardObserver_SpectralCurves(3,:) .* OFT_CurrentIlluminant_Spectrum_1nm_CIE31Range(2,:));
l_ZusedIllum = trapz(l_CIEStandardObserver_SpectralCurves(4,:) .* OFT_CurrentIlluminant_Spectrum_1nm_CIE31Range(2,:));

l_CCT = i_xy2cct([ l_XusedIllum / (l_XusedIllum + l_YusedIllum + l_ZusedIllum), l_YusedIllum / (l_XusedIllum + l_YusedIllum + l_ZusedIllum)]);

l_CCTIlluminant = HDM_OFT_GetBlackBodyRadiatorIllumination(l_CCT);

HDM_OFT_UI_PlotAndSave([OFT_CurrentIlluminant_Spectrum_1nm_CIE31Range(1,:); OFT_CurrentIlluminant_Spectrum_1nm_CIE31Range(2,:); ...
    l_CCTIlluminant(2, :); OFT_ReferenceIlluminant_Spectrum_1nm_CIE31Range(2, :)], ...
        'Normalized Scene Illumination Power Distribution', 'Wavelength (nm)', 'Normalized Radiance', ...
        {'Used Illuminant','CCT SPD for Used Illuminant','Target Illuminant'});


% parpool(2) disabled due to global env condition
% spmd

%% 4.7.2-4.7.4 reference tristimuli
HDM_OFT_Utils.OFT_DispSubTitle('4.7.2 - 4.7.4 prepare reference tristumuli');

%%
    
if (size(i_ReferenceObserver, 1) == 4 || isempty(strfind(i_ReferenceObserver, '.')) || not(isempty(strfind(i_ReferenceObserver, '.csv'))))

    [parOFT_PatchSetTristimuli_NeutralsCompensated, l_referenceWhite] = ...
        HDM_OFT_ComputeReferenceTristimuli4PatchSet(i_ReferenceObserver, OFT_ReferenceIlluminant_Spectrum_1nm_CIE31Range, OFT_PatchSet_SpectralCurve,...
        i_NeutralsCompensation, par4gOFT_w);

elseif(not(isempty(strfind(lower(i_CurrentObserver), '.tif'))) || not(isempty(strfind(lower(i_CurrentObserver), '.dpx'))))

    %//!!! not tested
    [parOFT_PatchSetTristimuli_NeutralsCompensated, o_IDT_b] = ...
        HDM_OFT_GetObservedTristimuliFromImage(i_ReferenceObserver, i_PreLinearisationCurve4CurrentObserverImage, i_indexOfReferenceWhite);            

    % todo reference white
    
end

%% 4.7.5 and 6 current tristumuli
HDM_OFT_Utils.OFT_DispSubTitle('4.7.5 and 6 prepare camera tristumuli');

if (size(i_CurrentObserver, 1) == 4 || isempty(strfind(i_CurrentObserver, '.')) || not(isempty(strfind(i_CurrentObserver, '.csv'))))

    [parOFT_PatchSetCameraTristimuli, o_IDT_b] = ...
        HDM_OFT_ComputeCurrentTristimuli4PatchSet(i_CurrentObserver, OFT_CurrentIlluminant_Spectrum_1nm_CIE31Range, OFT_ReferenceIlluminant_Spectrum_1nm_CIE31Range, OFT_PatchSet_SpectralCurve);

elseif(not(isempty(strfind(lower(i_CurrentObserver), '.tif'))) || not(isempty(strfind(lower(i_CurrentObserver), '.dpx'))))

    [parOFT_PatchSetCameraTristimuli, o_IDT_b] = ...
        HDM_OFT_GetObservedTristimuliFromImage(i_CurrentObserver, i_PreLinearisationCurve4CurrentObserverImage, i_indexOfReferenceWhite);            

end
% end
% delete(gcp)

%% 4.7.7 B estimation precise

HDM_OFT_Utils.OFT_DispSubTitle('4.7.7 B estimation');

if not(exist('i_ObjectReflectancesWeightings', 'var'))
    
    i_ObjectReflectancesWeightings = [];
    
end

% for constrained methods tests, if Color Checker is used
%i_IndiceesOfObjectReflectances4CorrectMappedColour = 19;

if not(exist('i_IndiceesOfObjectReflectances4CorrectMappedColour', 'var')) || isempty(i_IndiceesOfObjectReflectances4CorrectMappedColour)

    % contains tested ACES based method if i_ErrorMinimizationDomain is Lab
    %old script colorimetric transform matrix, i.e. left side before colour
    [o_IDT_B, l_resnormBEstimation] = ...
        HDM_OFT_Estimate3x3TransformationMatrixUnconstrained(i_ErrorMinimizationDomain, ...
            parOFT_PatchSetTristimuli_NeutralsCompensated', parOFT_PatchSetCameraTristimuli', l_M, l_referenceWhite, ...
            i_ObjectReflectancesWeightings);
    
else % constrained
    
    %new script colorimetric transform matrix, i.e. right side before colour^T
    [o_IDT_B, l_resnormBEstimation] = ...
        HDM_OFT_Estimate3x3TransformationMatrixConstrained(i_ErrorMinimizationDomain, ...
            parOFT_PatchSetTristimuli_NeutralsCompensated', parOFT_PatchSetCameraTristimuli', l_M, ...
            l_referenceWhite, i_ObjectReflectancesWeightings, i_IndiceesOfObjectReflectances4CorrectMappedColour);
       
    %so we must transpose for backward compat.
    o_IDT_B = o_IDT_B';
    
end

o_referenceWhite = l_referenceWhite;
o_resnormBEstimation = l_resnormBEstimation;

%% evaluate

l_transformedCurrentValues =  l_M * o_IDT_B * parOFT_PatchSetCameraTristimuli;
l_referenceValues = parOFT_PatchSetTristimuli_NeutralsCompensated;   

l_transformedCurrentValuesChromaticities = l_transformedCurrentValues ./ repmat(sum(l_transformedCurrentValues, 1), 3, 1);
l_referenceValuesChromaticities = l_referenceValues ./ repmat(sum(l_referenceValues, 1), 3, 1);

[OFT_CIE31_colorValuePartsWaveLength,...
OFT_CIE31_colorValueParts_x,OFT_CIE31_colorValueParts_y,OFT_CIE31_colorValueParts_z]=...
    HDM_OFT_CIEStandard.ColorValuePartsForSpectralColoursCurve(HDM_OFT_CIEStandard.StandardObserver1931_2Degrees);  

OFT_ColorCheckerPatchSetReference_xyY=HDM_OFT_PatchSet.GetCIE31_2Degress_D50_ColorChecker_BabelColorReferences();

l_plotArgs = {{[OFT_CIE31_colorValueParts_x; OFT_CIE31_colorValueParts_x(1)], [OFT_CIE31_colorValueParts_y;OFT_CIE31_colorValueParts_y(1)],'-'}; ...
              {l_referenceValuesChromaticities(1, :), l_referenceValuesChromaticities(2, :),'r+'}; ...  
              {l_transformedCurrentValuesChromaticities(1, :), l_transformedCurrentValuesChromaticities(2, :), 'bx'}; ...
              {l_referenceWhite(1, 1) / sum(l_referenceWhite(:)), l_referenceWhite(2, 1) / sum(l_referenceWhite(:)), 'gx'}};
          
l_labels2 = cellstr( num2str([1 : size(l_referenceValuesChromaticities, 2)]') );
l_dist4Label = 0.003;

l_textArgs = ...
    {{l_referenceValuesChromaticities(1, :) + l_dist4Label, l_referenceValuesChromaticities(2, :) + l_dist4Label, ...
    l_labels2, ...
    'VerticalAlignment','bottom', 'HorizontalAlignment','left'}; ...
    ...
    {l_transformedCurrentValuesChromaticities(1, :) - l_dist4Label, l_transformedCurrentValuesChromaticities(2, :) - l_dist4Label, ...
    l_labels2, ...
    'VerticalAlignment','top', 'HorizontalAlignment','right'}};

HDM_OFT_UI_PlotAndSave(l_plotArgs, ...
        'CIE 1931 2 Degree Standard Observer Chromaticities', 'x', 'y', ...
        {'spectral locus', 'reference', 'transformed', 'white'}, l_textArgs);

o_meanDeltaE2000 = -1;
o_stdDevDeltaE2000 = -1;
 
l_refLab = HDM_OFT_ColorConversions.OFT_CIELab(l_referenceValues, l_referenceWhite);
l_transformedCurLab = HDM_OFT_ColorConversions.OFT_CIELab(l_transformedCurrentValues, l_referenceWhite);

%% compute delta E 

l_deltaE2000 = deltaE2000(l_refLab', l_transformedCurLab'); 

l_cnames = {'L ref','a ref','b ref','L','a','b','delta E 2000'};
l_patchIndex = 1 : size(l_refLab, 2);
l_patchIndex = num2cell(l_patchIndex, size(l_refLab, 2));
HDM_OFT_UI_PlotTableAndSave...
    ([l_refLab', l_transformedCurLab', l_deltaE2000'], 'Training Set Ref vs. Transformed', l_cnames, l_patchIndex);

o_meanDeltaE2000 = mean(l_deltaE2000);
o_stdDevDeltaE2000 = std(l_deltaE2000);  

%%white check

[OFT_MRef,OFT_wRef]=HDM_OFT_IDT_ReferenceCamera.GetDefinition(i_ReferenceDomain);
OFT_Reference2StandardPrimaries=OFT_MRef;
OFT_MOverall=OFT_Reference2StandardPrimaries*o_IDT_B;

l_wU = o_IDT_B * [1; 1; 1];
l_wU = OFT_MOverall * [1; 1; 1];

whitePatchCameraRelatedE=parOFT_PatchSetCameraTristimuli(:,19);
whitePatchCameraRelatedE_bScaled=(o_IDT_b./min(o_IDT_b)).*whitePatchCameraRelatedE;
xyzNormScaled=whitePatchCameraRelatedE_bScaled(1)+whitePatchCameraRelatedE_bScaled(2)+whitePatchCameraRelatedE_bScaled(3);

wPxS=whitePatchCameraRelatedE_bScaled(1)/xyzNormScaled;
wPyS=whitePatchCameraRelatedE_bScaled(2)/xyzNormScaled;
wPzS=whitePatchCameraRelatedE_bScaled(3)/xyzNormScaled;

whitePatchCameraRelatedE_BConverted=OFT_MOverall*whitePatchCameraRelatedE_bScaled;
xyzNorm=whitePatchCameraRelatedE_BConverted(1)+whitePatchCameraRelatedE_BConverted(2)+whitePatchCameraRelatedE_BConverted(3);

wPx=whitePatchCameraRelatedE_BConverted(1)/xyzNorm;
wPy=whitePatchCameraRelatedE_BConverted(2)/xyzNorm;
wPz=whitePatchCameraRelatedE_BConverted(3)/xyzNorm;

HDM_OFT_Utils.OFT_DispTitle('idt profile successfully created');

end

function [o_PatchSetReferenceTristimuli, o_referenceWhite] = HDM_OFT_ComputeReferenceTristimuli4PatchSet...
    (i_ReferenceObserver, ...
    i_ReferenceIlluminant_Spectrum_1nm_CIE31Range, ...
    i_ObjectReflectances_SpectralCurve,...
    OFT_NeutralsCompensation, OFT_w)

%% CIE31 curves
HDM_OFT_Utils.OFT_DispSubTitle('setup CIE standard observers curves');

if size(i_ReferenceObserver, 1) == 4
    
    OFT_CIEStandardObserver_SpectralCurves = i_ReferenceObserver;
    
else
    
    OFT_CIEStandardObserver_SpectralCurves = HDM_OFT_CIEStandard.GetStandardObserverCurves(i_ReferenceObserver);

end

HDM_OFT_UI_PlotAndSave([OFT_CIEStandardObserver_SpectralCurves(1,:); ...
    OFT_CIEStandardObserver_SpectralCurves(2:4,:)], ...
        'Normalized Response of Standard Observer', 'Wavelength (nm)', 'Normalized Response', ...
        {'x','y','z'});

%% 4.7.2 patches
HDM_OFT_Utils.OFT_DispSubTitle('4.7.2 compute tristimuli for patches');
[OFT_PatchSetTristimuli, OFT_PatchSetTristimuli_Chromaticities]=...
        HDM_OFT_TristimuliCreator.CreateFromSpectrum(...
                OFT_CIEStandardObserver_SpectralCurves,...
                i_ReferenceIlluminant_Spectrum_1nm_CIE31Range,...
                i_ObjectReflectances_SpectralCurve);

%% plausibility check xyY //!!! for other patch sets must be ignored
HDM_OFT_Utils.OFT_DispSubTitle('xyY plausibility check for CIE1931 2 degrees With Babelcolor ColorChecker for D50');

[OFT_CIE31_colorValuePartsWaveLength,...
OFT_CIE31_colorValueParts_x,OFT_CIE31_colorValueParts_y,OFT_CIE31_colorValueParts_z]=...
HDM_OFT_CIEStandard.ColorValuePartsForSpectralColoursCurve(HDM_OFT_CIEStandard.StandardObserver1931_2Degrees);  
OFT_ColorCheckerPatchSetReference_xyY=HDM_OFT_PatchSet.GetCIE31_2Degress_D50_ColorChecker_BabelColorReferences();

l_plotArgs = {{[OFT_CIE31_colorValueParts_x;OFT_CIE31_colorValueParts_x(1)],[OFT_CIE31_colorValueParts_y;OFT_CIE31_colorValueParts_y(1)],'-'}; ...
              {OFT_PatchSetTristimuli_Chromaticities(1,:),OFT_PatchSetTristimuli_Chromaticities(2,:),'r+'}; ...  
              {OFT_ColorCheckerPatchSetReference_xyY(:,1), OFT_ColorCheckerPatchSetReference_xyY(:,2), 'bx'}};

HDM_OFT_UI_PlotAndSave(l_plotArgs, ...
        'CIE 1931 2 Degree Standard Observer Chromaticities', 'x', 'y');

%% 4.7.3 scene adopted white tristimulus, here the illumination source
HDM_OFT_Utils.OFT_DispSubTitle('4.7.3 setup tristimuli for scene adopetd white currently daylight from above used');
%figure

HDM_OFT_UI_PlotAndSave([i_ReferenceIlluminant_Spectrum_1nm_CIE31Range(1,:); ...
    i_ReferenceIlluminant_Spectrum_1nm_CIE31Range(2,:)], ...
        'Illuminant Spectrum', 'Wavelength (nm)', 'Normalized Spectral Power Distribution');

OFT_Illumination_Scale=1;
OFT_Illumination_Norm=1;

OFT_Xw = OFT_Illumination_Scale*trapz(OFT_CIEStandardObserver_SpectralCurves(2,:) .* i_ReferenceIlluminant_Spectrum_1nm_CIE31Range(2,:))/OFT_Illumination_Norm;
OFT_Yw = OFT_Illumination_Scale*trapz(OFT_CIEStandardObserver_SpectralCurves(3,:) .* i_ReferenceIlluminant_Spectrum_1nm_CIE31Range(2,:))/OFT_Illumination_Norm;
OFT_Zw = OFT_Illumination_Scale*trapz(OFT_CIEStandardObserver_SpectralCurves(4,:) .* i_ReferenceIlluminant_Spectrum_1nm_CIE31Range(2,:))/OFT_Illumination_Norm;

OFT_WwUnscaled=[OFT_Xw,OFT_Yw,OFT_Zw]';
OFT_Ww=OFT_WwUnscaled ./ OFT_Yw;

HDM_OFT_Utils.OFT_DispTitle('White XYZ plausibility check');
disp(OFT_Ww);
OFT_WwxyY=[OFT_Xw/(OFT_Xw + OFT_Yw + OFT_Zw),OFT_Yw/(OFT_Xw + OFT_Yw + OFT_Zw),OFT_Zw/(OFT_Xw + OFT_Yw + OFT_Zw)]';
disp(OFT_WwxyY);

o_referenceWhite = OFT_w;

%% 4.7.4 adjust tristimuli of training colours to compensate scene adopted
HDM_OFT_Utils.OFT_DispSubTitle('4.7.4 adjust tristimuli of training colours to compensate scene adopted');
o_PatchSetReferenceTristimuli=...
    HDM_OFT_ColorNeutralCompensations.OFT_CompensateTristimuliForDifferentWhite(OFT_NeutralsCompensation, OFT_PatchSetTristimuli, OFT_Ww, OFT_w);

end

function [OFT_PatchSetCameraTristimuli, OFT_IDT_b] = HDM_OFT_ComputeCurrentTristimuli4PatchSet...
    (OFT_MeasuredCameraResponse, ...
    OFT_Illuminant_Spectrum_1nm_CIE31Range, ...
    OFT_TargetIlluminant_Spectrum_1nm_CIE31Range, ...
    OFT_PatchSet_SpectralCurve)

    %not aequidistant
    
    HDM_OFT_UI_PlotAndSave([OFT_Illuminant_Spectrum_1nm_CIE31Range(1,:); ...
        OFT_Illuminant_Spectrum_1nm_CIE31Range(2,:)], ...
        'Normalized Scene Illumination Power Distribution', 'Wavelength (nm)', 'Normalized Spectral Power Distribution', ...
        {});
    
    HDM_OFT_Utils.OFT_DispSubTitle('start camera spectral response based tristimuli computation');    

    if size(OFT_MeasuredCameraResponse, 1) == 4

        OFT_CameraSpectralResponse_1nm_CIE31Range = OFT_MeasuredCameraResponse;

    else

        OFT_CameraSpectralResponse_1nm_CIE31Range = HDM_OFT_GetSpectralResponse(OFT_MeasuredCameraResponse);

    end;

    HDM_OFT_UI_PlotAndSave([OFT_CameraSpectralResponse_1nm_CIE31Range(1,:); ...
        OFT_CameraSpectralResponse_1nm_CIE31Range(4,:); ...
        OFT_CameraSpectralResponse_1nm_CIE31Range(3,:); ...
        OFT_CameraSpectralResponse_1nm_CIE31Range(2,:)], ...
        'Current Observer Spectral Response', 'Wavelength (nm)', 'Normalized Spectral Response', ...
        {'b','g','r'});

    %% 4.7.5 camera system white balance factors
    HDM_OFT_Utils.OFT_DispTitle('4.7.5 camera system white balance factors');
    
    HDM_OFT_Utils.OFT_DispTitle('4.7.6 compute white balanced linearized camera system response values of training colours');
    [OFT_PatchSetCameraTristimuli,OFT_PatchSetCameraTristimuli_ColorValueParts]=...
            HDM_OFT_TristimuliCreator.CreateFromSpectrum(...
                    OFT_CameraSpectralResponse_1nm_CIE31Range,...
                    ...
                    ... % ACES like: 
                    ... %HDM_OFT_GetIlluminantSpectrum('E'), ...
                    ...
                    ... %OFT modification
                    OFT_Illuminant_Spectrum_1nm_CIE31Range,... 
                    ...
                    OFT_PatchSet_SpectralCurve);         

    % ACES like: 
    OFT_CAM_Xw = trapz(OFT_CameraSpectralResponse_1nm_CIE31Range(2,:) .* OFT_TargetIlluminant_Spectrum_1nm_CIE31Range(2,:));
    OFT_CAM_Yw = trapz(OFT_CameraSpectralResponse_1nm_CIE31Range(3,:) .* OFT_TargetIlluminant_Spectrum_1nm_CIE31Range(2,:));
    OFT_CAM_Zw = trapz(OFT_CameraSpectralResponse_1nm_CIE31Range(4,:) .* OFT_TargetIlluminant_Spectrum_1nm_CIE31Range(2,:));

    % OFT modification
	OFT_CAM_Xw=trapz(OFT_CameraSpectralResponse_1nm_CIE31Range(2,:) .* OFT_Illuminant_Spectrum_1nm_CIE31Range(2,:));
	OFT_CAM_Yw=trapz(OFT_CameraSpectralResponse_1nm_CIE31Range(3,:) .* OFT_Illuminant_Spectrum_1nm_CIE31Range(2,:));
	OFT_CAM_Zw=trapz(OFT_CameraSpectralResponse_1nm_CIE31Range(4,:) .* OFT_Illuminant_Spectrum_1nm_CIE31Range(2,:));

    OFT_CAM_WwUnscaled = [OFT_CAM_Xw;OFT_CAM_Yw;OFT_CAM_Zw];               
    OFT_IDT_b = 1./OFT_CAM_WwUnscaled;

    %% 4.7.6 compute white balanced linearized camera system response values of training colours
    
    OFT_PatchSetCameraTristimuli3 = OFT_PatchSetCameraTristimuli;
    OFT_PatchSetCameraTristimuliW = OFT_PatchSetCameraTristimuli3 .* repmat(OFT_IDT_b,[1,size(OFT_PatchSetCameraTristimuli3,2)]);
    
    OFT_PatchSetCameraTristimuli=OFT_PatchSetCameraTristimuliW;

end

function [OFT_PatchSetCameraTristimuli, OFT_IDT_b] = HDM_OFT_GetObservedTristimuliFromImage...
    (i_Image, ...
    i_PreLinearisationCurve, ...
    i_whitePosition)

    %% annex B estimation by test chart image from camera
    HDM_OFT_Utils.OFT_DispTitle('start camera rgb image based tristimuli computation (Annex B)');    

    HDM_OFT_Utils.OFT_DispSubTitle('search for test chart in image');
    OFT_cameraImageOfTestChartOrigin=HDM_OFT_ImageExportImport.ImportImage(i_Image, i_PreLinearisationCurve);

    if usejava('Desktop')
        imshow(OFT_cameraImageOfTestChartOrigin);
    end

    OFT_cameraImageOfTestChart = double(OFT_cameraImageOfTestChartOrigin);
    [OFT_cameraImageOfTestChart_PatchLocations,OFT_cameraImageOfTestChart_PatchColours] = CCFind(OFT_cameraImageOfTestChart);
    OFT_PatchSetCameraTristimuli = OFT_cameraImageOfTestChart_PatchColours;

    %% 4.7.5 camera system white balance factors
    OFT_CAM_Ww = OFT_PatchSetCameraTristimuli(:, i_whitePosition);
    OFT_IDT_b=[1;1;1];%1./OFT_CAM_Ww;

    OFT_IDT_b=1./OFT_CAM_Ww;

    %% 4.7.6 compute white balanced linearized camera system response values of training colours

    OFT_PatchSetCameraTristimuli3 = OFT_PatchSetCameraTristimuli;
    OFT_PatchSetCameraTristimuliW = OFT_PatchSetCameraTristimuli3 .* repmat(OFT_IDT_b,[1,size(OFT_PatchSetCameraTristimuli3,2)]);

    OFT_PatchSetCameraTristimuli=OFT_PatchSetCameraTristimuliW;

end