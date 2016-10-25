function [OFT_IDT_File, OFT_IDT_B, OFT_IDT_b]=HDM_OFT_IDT_MinimumError...
    (OFT_In_PatchMeasurementFile, ...
    OFT_In_CameraMeasurementFile, ...
    OFT_In_PreLinearisationCurve, ...
    OFT_In_NeutralsCompensation, ...
    OFT_In_ErrorMinimizationDomain, ...
    OFT_In_IlluminantSpectrum,...
    OFT_In_UsedIlluminantSpectrum,...
    OFT_In_ReferenceDomain,...
    OFT_In_StandardObserver)

%% defaults
OFT_Env=HDM_OFT_InitEnvironment(); 

HDM_OFT_Utils.OFT_DispTitle('start IDT creation');

if(exist('OFT_In_PatchMeasurementFile','var')==0)
    disp('using reference patch mesurements');
    OFT_In_PatchMeasurementFile=HDM_OFT_PatchSet.GretagMacbethColorChecker();
else
    disp(OFT_In_PatchMeasurementFile);
end

if(exist('OFT_In_CameraMeasurementFile','var')==0)
    disp('using reference camera mesurements');
    OFT_MeasuredCameraResponseFileName=strcat(OFT_Env.OFT_RootDataDir,'/cameraMeasurementReference/arri_d21_spectral_response_02.csv');
    %OFT_MeasuredCameraResponseFileName=strcat(OFT_Env.OFT_RootDataDir,'/cameraImagesReference/sTake005_Img0000005.TIF');
else
    OFT_MeasuredCameraResponseFileName=OFT_In_CameraMeasurementFile;
end

if(exist('OFT_In_PreLinearisationCurve','var')==0)
    disp('using no linearization');
    OFT_In_PreLinearisationCurve='';
end

if(exist('OFT_In_NeutralsCompensation','var')==0)
    disp('using default neutral compensation');
    OFT_NeutralsCompensation=HDM_OFT_ColorNeutralCompensations.ChromaticAdaptationBradfordType();  
else
    OFT_NeutralsCompensation=OFT_In_NeutralsCompensation;  
end

if(exist('OFT_In_ErrorMinimizationDomain','var')==0)
    disp('using default neutral compensation');
    OFT_ErrorMinimizationDomain='XYZ';
else
    OFT_ErrorMinimizationDomain=OFT_In_ErrorMinimizationDomain;
    disp(OFT_In_ErrorMinimizationDomain);  
end

if(exist('OFT_In_IlluminantSpectrum','var')==0 || isempty(OFT_In_IlluminantSpectrum))
    disp('using default neutral compensation');   
    OFT_IlluminantSpectrum='D55';    
else
    OFT_IlluminantSpectrum=OFT_In_IlluminantSpectrum;  
end

if(exist('OFT_In_UsedIlluminantSpectrum','var')==0 || isempty(OFT_In_UsedIlluminantSpectrum))
    disp('using default neutral compensation');   
    OFT_UsedIlluminantSpectrum='D55';    
else
    OFT_UsedIlluminantSpectrum=OFT_In_UsedIlluminantSpectrum;  
end

if(exist('OFT_In_ReferenceDomain','var')==0)
    disp('using default reference domain');   
    OFT_ReferenceDomain=HDM_OFT_IDT_ReferenceCamera.CIEType();   
else
    OFT_ReferenceDomain=OFT_In_ReferenceDomain;  
end

if(exist('OFT_In_StandardObserver','var')==0)
    disp('using default standard');   
    OFT_StandardObserver=HDM_OFT_CIEStandard.StandardObserver1931_2Degrees();   
else
    OFT_StandardObserver=OFT_In_StandardObserver;  
end

%%reference camera
HDM_OFT_Utils.OFT_DispSubTitle('setup reference camera');

global gOFT_PatchSetTristimuli_NeutralsCompensated
global gOFT_PatchSetCameraTristimuli

%XYZ to RICD 
global gOFT_M;
global gOFT_w;
[gOFT_M,gOFT_w]=HDM_OFT_IDT_ReferenceCamera.GetDefinition(OFT_ReferenceDomain);

%% illuminat spectrum aquisition
HDM_OFT_Utils.OFT_DispSubTitle('illuminat spectrum aquisition');
OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetIlluminantSpectrum(OFT_IlluminantSpectrum);
OFT_UsedIlluminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetIlluminantSpectrum(OFT_UsedIlluminantSpectrum);

%%patches spectrum aquisition
%10 nm resolution
HDM_OFT_Utils.OFT_DispSubTitle('read patch spectra');
OFT_PatchSet_SpectralCurve=HDM_OFT_PatchSet.GetPatchSpectra(OFT_In_PatchMeasurementFile);

par4gOFT_w=gOFT_w;

%% visualize spectral refelctance subset for colorchecker

figure ('Name','patch reflectances for selected improved colours')
plot(OFT_PatchSet_SpectralCurve(1,:),OFT_PatchSet_SpectralCurve(2,:),...
OFT_PatchSet_SpectralCurve(1,:),OFT_PatchSet_SpectralCurve(7,:),...
OFT_PatchSet_SpectralCurve(1,:),OFT_PatchSet_SpectralCurve(19,:))
%semilogy(OFT_CIEStandardObserver_SpectralCurves(1,:),OFT_CIEStandardObserver_SpectralCurves(2:4,:))
xlabel('Wavelength (nm)')
ylabel('Normalized Reflectance')
legend({'Dark Skin','Orange','Cyan'})
grid on
grid minor

%% visualize SPD of used illuminant, CCT related spectrum and target illuminant

l_CIEStandardObserver_SpectralCurves=HDM_OFT_CIEStandard.GetStandardObserverCurves(OFT_StandardObserver);

l_XusedIllum = trapz(l_CIEStandardObserver_SpectralCurves(2,:) .* OFT_UsedIlluminant_Spectrum_1nm_CIE31Range(2,:));
l_YusedIllum = trapz(l_CIEStandardObserver_SpectralCurves(3,:) .* OFT_UsedIlluminant_Spectrum_1nm_CIE31Range(2,:));
l_ZusedIllum = trapz(l_CIEStandardObserver_SpectralCurves(4,:) .* OFT_UsedIlluminant_Spectrum_1nm_CIE31Range(2,:));


l_CCT = i_xy2cct([ l_XusedIllum / (l_XusedIllum + l_YusedIllum + l_ZusedIllum), l_YusedIllum / (l_XusedIllum + l_YusedIllum + l_ZusedIllum)]);

l_CCTIlluminant = HDM_OFT_GetBlackBodyRadiatorIllumination(l_CCT);

figure ('Name','normalized scene illumination power distribution')
plot(OFT_UsedIlluminant_Spectrum_1nm_CIE31Range(1,:), OFT_UsedIlluminant_Spectrum_1nm_CIE31Range(2, :),...
l_CCTIlluminant(1, :), l_CCTIlluminant(2, :), ...
OFT_Illuminant_Spectrum_1nm_CIE31Range(1,:), OFT_Illuminant_Spectrum_1nm_CIE31Range(2, :))
xlabel('Wavelength (nm)')
ylabel('Normalized Radiance')
legend({'Used Illuminant','CCT SPD for Used Illuminant','Target Illuminant'})
grid on
grid minor


% parpool(2) disabled due to global env condition
% spmd

%% 4.7.2-4.7.4 reference tristimuli
HDM_OFT_Utils.OFT_DispSubTitle('4.7.2 - 4.7.4 prepare reference tristumuli');
[parOFT_PatchSetTristimuli_NeutralsCompensated,referenceWhite] = ComputeReferenceTristimuli4PatchSet...
    (OFT_StandardObserver, OFT_Illuminant_Spectrum_1nm_CIE31Range,OFT_PatchSet_SpectralCurve,...%HDM_OFT_GetIlluminantSpectrum('D50') ,OFT_PatchSet_SpectralCurve,...//!!! D50 for fix for test cause colorchecker Lab we have only for D50
    OFT_NeutralsCompensation, par4gOFT_w);

%% 4.7.5 and 6 camera tristumuli
HDM_OFT_Utils.OFT_DispSubTitle('4.7.5 and 6 prepare camera tristumuli');
[parOFT_PatchSetCameraTristimuli, OFT_IDT_b] = ComputeCameraTristimuli4PatchSet...
    (OFT_MeasuredCameraResponseFileName, OFT_In_PreLinearisationCurve, OFT_UsedIlluminant_Spectrum_1nm_CIE31Range, OFT_Illuminant_Spectrum_1nm_CIE31Range, OFT_PatchSet_SpectralCurve);

% end
% delete(gcp)

gOFT_PatchSetTristimuli_NeutralsCompensated=parOFT_PatchSetTristimuli_NeutralsCompensated;

if(strcmp(OFT_NeutralsCompensation, HDM_OFT_ColorNeutralCompensations.NoneType()))
    gOFT_w=referenceWhite;
end

gOFT_PatchSetCameraTristimuli=parOFT_PatchSetCameraTristimuli;

%% 4.7.7 B estimation precise
%//!!!todo weight implementation
HDM_OFT_Utils.OFT_DispSubTitle('4.7.7 B estimation');
[OFT_IDT_B, OFT_resnormBEstimation]=EstimateIDTMatrix(OFT_ErrorMinimizationDomain);

disp('reference matrix');
gOFT_M

disp('estimated matrix');
OFT_IDT_B

disp('reference matrix by estimated matrix');
gOFT_M*OFT_IDT_B

%% write idt profile file and append results to statistics
HDM_OFT_Utils.OFT_DispTitle('write idt profile file and append results to statistics');
OFT_IDT_File = WriteIDTProfileAndStatEntry...
    (OFT_Env, OFT_MeasuredCameraResponseFileName, OFT_resnormBEstimation, OFT_IDT_B, OFT_IDT_b,...
    OFT_NeutralsCompensation, OFT_IlluminantSpectrum, OFT_ErrorMinimizationDomain);

%%white check
[OFT_MRef,OFT_wRef]=HDM_OFT_IDT_ReferenceCamera.GetDefinition(OFT_ReferenceDomain);
OFT_Reference2StandardPrimaries=OFT_MRef;
OFT_MOverall=OFT_Reference2StandardPrimaries*OFT_IDT_B;

whitePatchCameraRelatedE=parOFT_PatchSetCameraTristimuli(:,19);

whitePatchCameraRelatedE_bScaled=(OFT_IDT_b./min(OFT_IDT_b)).*whitePatchCameraRelatedE;

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

function [OFT_PatchSetReferenceTristimuli, referenceWhite] = ComputeReferenceTristimuli4PatchSet...
    (OFT_StandardObserver, OFT_Illuminant_Spectrum_1nm_CIE31Range, OFT_PatchSet_SpectralCurve,...
    OFT_NeutralsCompensation, OFT_w)

%% CIE31 curves
HDM_OFT_Utils.OFT_DispSubTitle('setup CIE standard observers curves');
OFT_CIEStandardObserver_SpectralCurves=HDM_OFT_CIEStandard.GetStandardObserverCurves(OFT_StandardObserver);

figure 
subplot(2,2,1)
plot(OFT_CIEStandardObserver_SpectralCurves(1,:),OFT_CIEStandardObserver_SpectralCurves(2:4,:))
%semilogy(OFT_CIEStandardObserver_SpectralCurves(1,:),OFT_CIEStandardObserver_SpectralCurves(2:4,:))
xlabel('wavelength in nm')
ylabel('relative sensitivity of standard observer')
legend({'x','y','z'})
title(OFT_StandardObserver);

%% 4.7.2 patches
HDM_OFT_Utils.OFT_DispSubTitle('4.7.2 compute tristimuli for patches');
[OFT_PatchSetTristimuli,OFT_PatchSetTristimuli_ColorValueParts]=...
        HDM_OFT_TristimuliCreator.CreateFromSpectrum(...
                OFT_CIEStandardObserver_SpectralCurves,...
                OFT_Illuminant_Spectrum_1nm_CIE31Range,...
                OFT_PatchSet_SpectralCurve);

%% plausibility check xyY //!!! for other patch sets must be ignored
HDM_OFT_Utils.OFT_DispSubTitle('xyY plausibility check for CIE1931 2 degress With Babelcolor ColorChecker for D50');

[OFT_CIE31_colorValuePartsWaveLength,...
OFT_CIE31_colorValueParts_x,OFT_CIE31_colorValueParts_y,OFT_CIE31_colorValueParts_z]=...
HDM_OFT_CIEStandard.ColorValuePartsForSpectralColoursCurve(HDM_OFT_CIEStandard.StandardObserver1931_2Degrees);  
OFT_ColorCheckerPatchSetReference_xyY=HDM_OFT_PatchSet.GetCIE31_2Degress_D50_ColorChecker_BabelColorReferences();

subplot(2,2,2)
%figure
%scatter(OFT_PatchSetTristimuli_ColorValueParts(1,:),OFT_PatchSetTristimuli_ColorValueParts(2,:))
plot([OFT_CIE31_colorValueParts_x;OFT_CIE31_colorValueParts_x(1)],[OFT_CIE31_colorValueParts_y;OFT_CIE31_colorValueParts_y(1)],'-',...
    OFT_PatchSetTristimuli_ColorValueParts(1,:),OFT_PatchSetTristimuli_ColorValueParts(2,:),'r+',...
    OFT_ColorCheckerPatchSetReference_xyY(:,1),OFT_ColorCheckerPatchSetReference_xyY(:,2),'bx')

xlabel('x')
ylabel('y')
title('CIE31 x y color value parts');

%% 4.7.3 scene adopted white tristimulus, here the illumination source
HDM_OFT_Utils.OFT_DispSubTitle('4.7.3 setup tristimuli for scene adopetd white currently daylight from above used');
%figure

subplot(2,2,3)
plot(OFT_Illuminant_Spectrum_1nm_CIE31Range(1,:),OFT_Illuminant_Spectrum_1nm_CIE31Range(2,:))
xlabel('wavelength in nm')
ylabel('relative illumination power distribution')
title('Daylight for CIE31 range');

OFT_Illumination_Scale=1;
OFT_Illumination_Norm=1;

OFT_Xw=OFT_Illumination_Scale*trapz(OFT_CIEStandardObserver_SpectralCurves(2,:) .* OFT_Illuminant_Spectrum_1nm_CIE31Range(2,:))/OFT_Illumination_Norm;
OFT_Yw=OFT_Illumination_Scale*trapz(OFT_CIEStandardObserver_SpectralCurves(3,:) .* OFT_Illuminant_Spectrum_1nm_CIE31Range(2,:))/OFT_Illumination_Norm;
OFT_Zw=OFT_Illumination_Scale*trapz(OFT_CIEStandardObserver_SpectralCurves(4,:) .* OFT_Illuminant_Spectrum_1nm_CIE31Range(2,:))/OFT_Illumination_Norm;

OFT_WwUnscaled=[OFT_Xw,OFT_Yw,OFT_Zw]';
OFT_Ww=100*(OFT_WwUnscaled./OFT_WwUnscaled(2));
HDM_OFT_Utils.OFT_DispTitle('Daylight XYZ plausibility check');
disp(OFT_Ww);

OFT_WwxyY=[OFT_Xw/(OFT_Xw + OFT_Yw + OFT_Zw),OFT_Yw/(OFT_Xw + OFT_Yw + OFT_Zw),OFT_Zw/(OFT_Xw + OFT_Yw + OFT_Zw)]';
disp(OFT_WwxyY);


OFT_PatchSetTristimuliNorm=100*(OFT_PatchSetTristimuli./OFT_WwUnscaled(2));

OFT_PatchSetTristimuli=OFT_PatchSetTristimuliNorm;
referenceWhite=OFT_Ww;

%% 4.7.4 adjust tristimuli of training colours to compensate scene adopted
HDM_OFT_Utils.OFT_DispSubTitle('4.7.4 adjust tristimuli of training colours to compensate scene adopted');
OFT_PatchSetReferenceTristimuli=...
    HDM_OFT_ColorNeutralCompensations.OFT_CompensateTristimuliForDifferentWhite(OFT_NeutralsCompensation, OFT_PatchSetTristimuli, OFT_Ww, OFT_w);

end

function [OFT_PatchSetCameraTristimuli, OFT_IDT_b] = ComputeCameraTristimuli4PatchSet...
    (OFT_MeasuredCameraResponseFileName, OFT_PreLinearisationCurve, OFT_Illuminant_Spectrum_1nm_CIE31Range, OFT_TargetIlluminant_Spectrum_1nm_CIE31Range, OFT_PatchSet_SpectralCurve)

    %not aequidistant
    
    subplot(2,2,2)
    plot(OFT_Illuminant_Spectrum_1nm_CIE31Range(1,:),OFT_Illuminant_Spectrum_1nm_CIE31Range(2,:))
    xlabel('wavelength in nm')
    ylabel('relative power distribution')
    title('relative scene illumination power distribution');

    if(strfind(OFT_MeasuredCameraResponseFileName, '.csv'))

        HDM_OFT_Utils.OFT_DispSubTitle('start camera spectral response based tristimuli computation');    

        OFT_CameraSpectralResponse_1nm_CIE31Range = HDM_OFT_GetSpectralResponse(OFT_MeasuredCameraResponseFileName);

        subplot(2,2,4)
        plot(OFT_CameraSpectralResponse_1nm_CIE31Range(1,:),OFT_CameraSpectralResponse_1nm_CIE31Range(4,:),...
            OFT_CameraSpectralResponse_1nm_CIE31Range(1,:),OFT_CameraSpectralResponse_1nm_CIE31Range(3,:),...
            OFT_CameraSpectralResponse_1nm_CIE31Range(1,:),OFT_CameraSpectralResponse_1nm_CIE31Range(2,:));
        xlabel('wavelength in nm')
        ylabel('relative spectral response')
        legend({'b','g','r'})
        %title('spectral response of camera system');
        
        figure
        plot(OFT_CameraSpectralResponse_1nm_CIE31Range(1,:),OFT_CameraSpectralResponse_1nm_CIE31Range(4,:),...
            OFT_CameraSpectralResponse_1nm_CIE31Range(1,:),OFT_CameraSpectralResponse_1nm_CIE31Range(3,:),...
            OFT_CameraSpectralResponse_1nm_CIE31Range(1,:),OFT_CameraSpectralResponse_1nm_CIE31Range(2,:));
        xlabel('wavelength in nm')
        ylabel('relative spectral response')
        legend({'b','g','r'})
        grid on
        grid minor
        %title('spectral response of camera system');

        %% 4.7.5 camera system white balance factors
        HDM_OFT_Utils.OFT_DispTitle('4.7.5 camera system white balance factors');
        
        OFT_CAM_Xw=trapz(OFT_CameraSpectralResponse_1nm_CIE31Range(2,:) .* OFT_TargetIlluminant_Spectrum_1nm_CIE31Range(2,:));
        OFT_CAM_Yw=trapz(OFT_CameraSpectralResponse_1nm_CIE31Range(3,:) .* OFT_TargetIlluminant_Spectrum_1nm_CIE31Range(2,:));
        OFT_CAM_Zw=trapz(OFT_CameraSpectralResponse_1nm_CIE31Range(4,:) .* OFT_TargetIlluminant_Spectrum_1nm_CIE31Range(2,:));

        OFT_CAM_WwUnscaled = [OFT_CAM_Xw;OFT_CAM_Yw;OFT_CAM_Zw];               
        OFT_IDT_b=1./OFT_CAM_WwUnscaled;

        %% 4.7.6 compute white balanced linearized camera system response values of training colours
        HDM_OFT_Utils.OFT_DispTitle('4.7.6 compute white balanced linearized camera system response values of training colours');
        [OFT_PatchSetCameraTristimuli,OFT_PatchSetCameraTristimuli_ColorValueParts]=...
                HDM_OFT_TristimuliCreator.CreateFromSpectrum(...
                        OFT_CameraSpectralResponse_1nm_CIE31Range,...
                        OFT_Illuminant_Spectrum_1nm_CIE31Range,...
                        OFT_PatchSet_SpectralCurve);         

        OFT_PatchSetCameraTristimuli3 = OFT_PatchSetCameraTristimuli;
        OFT_PatchSetCameraTristimuliW = OFT_PatchSetCameraTristimuli3 .* repmat(OFT_IDT_b,[1,size(OFT_PatchSetCameraTristimuli3,2)]);
        
        OFT_PatchSetCameraTristimuli=OFT_PatchSetCameraTristimuliW;

    elseif(isempty(strfind(lower(OFT_MeasuredCameraResponseFileName), '.tif'))||isempty(strfind(lower(OFT_MeasuredCameraResponseFileName), '.dpx')))%//!!!

        %% annex B estimation by test chart image from camera
        HDM_OFT_Utils.OFT_DispTitle('start camera rgb image based tristimuli computation (Annex B)');    

        HDM_OFT_Utils.OFT_DispSubTitle('search for test chart in image');
        OFT_cameraImageOfTestChartOrigin=HDM_OFT_ImageExportImport.ImportImage(OFT_MeasuredCameraResponseFileName, OFT_PreLinearisationCurve);
        
        if usejava('Desktop')
        	subplot(2,2,4) ,imshow(OFT_cameraImageOfTestChartOrigin);
        end
        
        OFT_cameraImageOfTestChart = double(OFT_cameraImageOfTestChartOrigin);
        [OFT_cameraImageOfTestChart_PatchLocations,OFT_cameraImageOfTestChart_PatchColours] = CCFind(OFT_cameraImageOfTestChart);
        OFT_PatchSetCameraTristimuli = OFT_cameraImageOfTestChart_PatchColours;
        
        %% 4.7.5 camera system white balance factors
        OFT_CAM_Ww = OFT_PatchSetCameraTristimuli(:,19);
        OFT_IDT_b=[1;1;1];%1./OFT_CAM_Ww;
        
        OFT_IDT_b=1./OFT_CAM_Ww;
        
        %% 4.7.6 compute white balanced linearized camera system response values of training colours

        OFT_PatchSetCameraTristimuli3 = OFT_PatchSetCameraTristimuli;
        OFT_PatchSetCameraTristimuliW = OFT_PatchSetCameraTristimuli3 .* repmat(OFT_IDT_b,[1,size(OFT_PatchSetCameraTristimuli3,2)]);
        
        OFT_PatchSetCameraTristimuli=OFT_PatchSetCameraTristimuliW;
        

    end

end

function [OFT_IDT_B,OFT_resnormBEstimation]=EstimateIDTMatrix(OFT_In_ErrorMinimizationDomain)

    OFT_IDT_BStart = ...
        [1 0 0;
        0 1 0;
        0 0 1];

    switch OFT_In_ErrorMinimizationDomain
        case 'Lab'
            [OFT_IDT_B,OFT_resnormBEstimation] = lsqnonlin(@OFT_IDT_MeritFunctionCoreLab,OFT_IDT_BStart);
        case 'Luv'
            [OFT_IDT_B,OFT_resnormBEstimation] = lsqnonlin(@OFT_IDT_MeritFunctionCoreLuv,OFT_IDT_BStart);
        case 'XYZ'
            [OFT_IDT_B,OFT_resnormBEstimation] = lsqnonlin(@OFT_IDT_MeritFunctionCoreXYZ,OFT_IDT_BStart);
        otherwise
    end

end

function F = OFT_IDT_MeritFunctionCoreLuv(B0)

global gOFT_M
global gOFT_w
global gOFT_PatchSetTristimuli_NeutralsCompensated
global gOFT_PatchSetCameraTristimuli

k = 1:size(gOFT_PatchSetTristimuli_NeutralsCompensated,2);

F = HDM_OFT_ColorConversions.OFT_CIELuv(gOFT_PatchSetTristimuli_NeutralsCompensated(:,k),gOFT_w)-...
    HDM_OFT_ColorConversions.OFT_CIELuv(gOFT_M*B0*gOFT_PatchSetCameraTristimuli(:,k),gOFT_w);

end

function F = OFT_IDT_MeritFunctionCoreLab(B0)

global gOFT_M
global gOFT_w
global gOFT_PatchSetTristimuli_NeutralsCompensated
global gOFT_PatchSetCameraTristimuli

k = 1:size(gOFT_PatchSetTristimuli_NeutralsCompensated,2);

F = HDM_OFT_ColorConversions.OFT_CIELab(gOFT_PatchSetTristimuli_NeutralsCompensated(:,k),gOFT_w)-...
    HDM_OFT_ColorConversions.OFT_CIELab(gOFT_M*B0*gOFT_PatchSetCameraTristimuli(:,k),gOFT_w);

end

function F = OFT_IDT_MeritFunctionCoreXYZ(B0)

global gOFT_PatchSetTristimuli_NeutralsCompensated
global gOFT_PatchSetCameraTristimuli

l=0;
u=0;
k = 1+l:(size(gOFT_PatchSetTristimuli_NeutralsCompensated,2)-u);

F = gOFT_PatchSetTristimuli_NeutralsCompensated(:,k)-B0*gOFT_PatchSetCameraTristimuli(:,k);

end

function OFT_IDTProfileFileName = WriteIDTProfileAndStatEntry...
    (OFT_Env, OFT_MeasuredCameraResponseFileName, OFT_resnormBEstimation, OFT_IDT_B, OFT_IDT_b,...
    OFT_NeutralsCompensation, OFT_IlluminantSpectrum, OFT_ErrorMinimizationDomain)

    if(isempty(strfind(OFT_IlluminantSpectrum,'.')))
        IlluminantStr=strcat('Illuminant_',OFT_IlluminantSpectrum);
    else
        [OFT_IlluminantSpectrumPath,OFT_IlluminantSpectrumName,OFT_IlluminantSpectrumExt] = fileparts(OFT_IlluminantSpectrum);
        IlluminantStr=strcat('Illuminant_FromFile_',OFT_IlluminantSpectrumName);
    end

    fin = fopen(strcat(OFT_Env.OFT_ConstraintsPath,'/IDT_Template.txt'));
    idtCreationDateStr=datestr(now,'yyyy-mm-dd_HH.MM.SS.FFF');
    OFT_IDT_File=strcat(OFT_Env.OFT_ProcessPath,'/IDT_',IlluminantStr,'_',idtCreationDateStr,'.ctl');
    fout = fopen(OFT_IDT_File,'wt');

    if ~exist(strcat(OFT_Env.OFT_StatisticsPath,'/IDTStat.csv'),'file')
        foutStat = fopen(strcat(OFT_Env.OFT_StatisticsPath,'/IDTStat.csv'),'wt');
        fprintf(foutStat,'idt file\t\t , measurement file\t\t , resnorm , B11 , B12 , B13 , B21 , B22 , B23 , B31 , B32 , B33 , scene adopted white, neutrals compensation, colour domain for error minimization\n');
    else
        foutStat = fopen(strcat(OFT_Env.OFT_StatisticsPath,'/IDTStat.csv'),'at');
    end

    fprintf(foutStat,'%s\t , ', strcat('IDT_',IlluminantStr,'_',idtCreationDateStr,'.ctl'));
    [~,OFT_MeasuredCameraResponseName,OFT_MeasuredCameraResponseExt] = fileparts(OFT_MeasuredCameraResponseFileName);
    fprintf(foutStat,'%s\t , ', strcat(OFT_MeasuredCameraResponseName,OFT_MeasuredCameraResponseExt));
    fprintf(foutStat,'%e , ', OFT_resnormBEstimation);

    while ~feof(fin)
       S = fgetl(fin);
       %s = strrep(s, '118520', '118521');

       if(strfind(S, '%'))

            if(strfind(S, '%IDT_DATE%'))

               fprintf(fout,'// Creation Date: %s\n',idtCreationDateStr);

            elseif (strfind(S, '%IDT_ILLUMINANT%'))

               fprintf(fout,'// Illuminant %s\n', OFT_IlluminantSpectrum);

            elseif(strfind(S, '%const float b['))

               fprintf(fout,'\tconst float b[] = { %f, %f, %f };\n',OFT_IDT_b(1),OFT_IDT_b(2),OFT_IDT_b(3));

            elseif(strfind(S, '%const float B1'))

               fprintf(fout,'\tconst float B[][] =     { { %f, %f, %f },\n',OFT_IDT_B(1,1),OFT_IDT_B(1,2),OFT_IDT_B(1,3));

               fprintf(foutStat,'%f , %f , %f , ',OFT_IDT_B(1,1),OFT_IDT_B(1,2),OFT_IDT_B(1,3));

            elseif(strfind(S, '%const float B2'))

               fprintf(fout,'\t\t\t  { %f, %f, %f },\n',OFT_IDT_B(2,1),OFT_IDT_B(2,2),OFT_IDT_B(2,3));

               fprintf(foutStat,'%f , %f , %f , ',OFT_IDT_B(2,1),OFT_IDT_B(2,2),OFT_IDT_B(2,3));

            elseif(strfind(S, '%const float B3'))

               fprintf(fout,'\t\t\t  { %f, %f, %f }};\n',OFT_IDT_B(3,1),OFT_IDT_B(3,2),OFT_IDT_B(3,3));

               fprintf(foutStat,'%f , %f , %f , ',OFT_IDT_B(3,1),OFT_IDT_B(3,2),OFT_IDT_B(3,3));

               fprintf(foutStat,'%s , %s , %s', OFT_IlluminantSpectrum, OFT_NeutralsCompensation, OFT_ErrorMinimizationDomain);

           end

       else

            fprintf(fout,'%s\n',S);

       end
    end

    fprintf(foutStat,'\n');

    fclose(fin);
    fclose(fout);
    fclose(foutStat);
    
    OFT_IDTProfileFileName=OFT_IDT_File;

end

% next two methods by
% https://de.mathworks.com/matlabcentral/fileexchange/13788-optprop-a-color-properties-toolbox/content/optprop/i_conversions/i_xy2cct.m?requestedDomain=www.mathworks.com

function z=i_xy2cct(xy)
%I_XYZ2CCT Calculate correlated color temperature.
%   CCT=I_XY2CCT(XYZ) with size(XYZ)=[M N ... P 2] returns
%   matrix CCT with size [M N ... P].
%
%   Remark:
%      This is a low level function, that has a rigid parameter passing
%      mechanism and no error handling. It is only to be used when the
%      need for speed is imperative. In all other cases,
%      use XY2CCT instead.
%
%	Example:
%      Find the correlated color temperatur for D65 with 2 degrees observer:
%
%         i_xy2cct(i_xyz2xy(wpt('D65/2'))) % returns CCT for D65
%
%   See: ROO2CCT, BLACKBODY, DILL

%   Wrapper and transcriber: Jerker Wgberg, 2005-03-30
%   Author: Neil Okamoto
%   Original at http://www.efg2.com/Lab/Library/UseNet/2001/0714.txt
%

% Part of the OptProp toolbox, $Version: 2.1 $
% Author:  Jerker Wgberg, More Research & DPC, Sweden
% Email: ['jerker.wagberg' char(64) 'more.se']

% $Id: i_xy2cct.m 23 2007-01-28 22:55:34Z jerkerw $

	% convert (x,y) to CIE 1960 (u,v)
	denom=-xy(:,1) + 6*xy(:,2) + 1.5;
	us = (2*xy(:,1)) ./ denom;
	vs = (3*xy(:,2)) ./ denom;

	z=cct(us,vs);
    
end

function t=cct(us,vs)
	% Computes correlated color temperature based on Robertson's method.
	% (cf. Wyszecki & Stiles, p.224-9)
	%
	persistent mirek ut vt tt
	if isempty(mirek)
		% ut(20) changed from W&S 0.24702 to 0.24792
		isodata=[ ...
			0,    0.18006,  0.26352,   -0.24341
		   10,    0.18066,  0.26589,   -0.25479
		   20,    0.18133,  0.26846,   -0.26876
		   30,    0.18208,  0.27119,   -0.28539
		   40,    0.18293,  0.27407,   -0.30470
		   50,    0.18388,  0.27709,   -0.32675
		   60,    0.18494,  0.28021,   -0.35156
		   70,    0.18611,  0.28342,   -0.37915
		   80,    0.18740,  0.28668,   -0.40955
		   90,    0.18880,  0.28997,   -0.44278
		   100,   0.19032,  0.29326,   -0.47888
		   125,   0.19462,  0.30141,   -0.58204
		   150,   0.19962,  0.30921,   -0.70471
		   175,   0.20525,  0.31647,   -0.84901
		   200,   0.21142,  0.32312,   -1.0182 
		   225,   0.21807,  0.32909,   -1.2168 
		   250,   0.22511,  0.33439,   -1.4512 
		   275,   0.23247,  0.33904,   -1.7298 
		   300,   0.24010,  0.34308,   -2.0637 
		   325,   0.24792,  0.34655,   -2.4681 
		   350,   0.25591,  0.34951,   -2.9641 
		   375,   0.26400,  0.35200,   -3.5814 
		   400,   0.27218,  0.35407,   -4.3633 
		   425,   0.28039,  0.35577,   -5.3762 
		   450,   0.28863,  0.35714,   -6.7262 
		   475,   0.29685,  0.35823,   -8.5955 
		   500,   0.30505,  0.35907,  -11.324 
		   525,   0.31320,  0.35968,  -15.628 
		   550,   0.32129,  0.36011,  -23.325 
		   575,   0.32931,  0.36038,  -40.770 
		   600,   0.33724,  0.36051, -116.45
		   ];
		mirek=isodata(:,1);	% temp (in microreciprocal kelvin)
		ut=isodata(:,2);	% u coord of intersection w/ blackbody locus
		vt=isodata(:,3);	% v coord of intersection w/ blackbody locus
		tt=isodata(:,4);	% slope of isotemp. line
		end

	% search for closest isotemp lines
	t=nan*zeros(size(us));
	ix=true(size(us));
	j=1;
	while any(ix) && j<=length(mirek)
		% dj = distance from (us,vs) to this isotemp line */
		dj(ix) = ((vs(ix) - vt(j)) - tt(j) * (us(ix) - ut(j))); %#ok<AGROW>

		ixe=false(size(ix));
		if j>1
			% we stop when di and dj changes sign, because this means we have
			% found isotemp lines that "straddle" our point.
			ixe(ix)=j>1 & (di(ix)<0 & dj(ix)>0 | di(ix)>0 & dj(ix)<0);
			dit=di(ixe)/sqrt(1+tt(j-1).^2);
			djt=dj(ixe)/sqrt(1+tt(j).^2);
			t(ixe)=1000000.0 ./ (mirek(j-1) + (dit ./ (dit - djt)) * (mirek(j) - mirek(j-1)));
			ix(ixe)=false;
			end
		di(ix) = dj(ix); %#ok<AGROW>
		j=j+1;
		end
end


