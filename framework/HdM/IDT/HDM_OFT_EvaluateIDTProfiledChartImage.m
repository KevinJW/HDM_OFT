function OFT_TransformedImage2View=HDM_OFT_EvaluateIDTProfiledChartImage...
    (OFT_IDT_B, OFT_IDT_b, ...
    OFT_In_CameraMeasurement, OFT_In_PreLinearisationCurve, ...
    OFT_In_IlluminationSpectrum, OFT_In_SceneIllumination, ... 
    OFT_In_PatchMeasurementFile, ...
    OFT_In_SpectralResponseFile, ...
    OFT_In_ReferenceDomain, OFT_In_NeutralsCompensation, OFT_In_StandardObserver, OFT_In_ReferenceImage2View)

OFT_Env=HDM_OFT_InitEnvironment();
HDM_OFT_Utils.OFT_DispTitle('evaluate transformed chart image');

if(exist('OFT_IDT_B','var')==0)
    
    close all;
    
    disp('using default matrix');
    OFT_B=  [1,  0,   0;...
            0,   1,   0;...
            0,   0,   1];
else   
    disp('using matrix');
    OFT_B=OFT_IDT_B
end

if(exist('OFT_IDT_b','var')==0)
    disp('using default white balance');
    OFT_b=  [1;...
             1;...
             1];
else   
    disp('using white balance');
    OFT_b=OFT_IDT_b
end

if(exist('OFT_In_CameraMeasurement','var')==0 || strcmp(OFT_In_CameraMeasurement,''))
    OFT_In_CameraMeasurement=strcat(OFT_Env.OFT_RootDataDir,'/cameraImagesReference/testimage-mmCenter.tif');% uint8 industry cam sample
    %OFT_In_CameraMeasurement=strcat(OFT_Env.OFT_RootDataDir,'/cameraImagesReference/Alexa_Wide_wb5000.tiff');% 32 bit float sample
    
    disp('using reference image'); 
    disp(OFT_In_CameraMeasurement);
else
    disp('using image'); 
    disp(OFT_In_CameraMeasurement);
end

if(exist('OFT_In_SpectralResponseFile','var')==0 || strcmp(OFT_In_SpectralResponseFile,''))
    %OFT_SpectralResponseFile = strcat(OFT_Env.OFT_RootDataDir,'/cameraResponsesReference/industryCAM.csv');
    OFT_SpectralResponseFile = strcat(OFT_Env.OFT_RootDataDir,'/cameraImagesReference/testimage-mmCenter.tif');
    
    disp('using reference camera response'); 
    disp(OFT_SpectralResponseFile);
else
    OFT_SpectralResponseFile=OFT_In_SpectralResponseFile;
    disp('using given camera response'); 
    disp(OFT_SpectralResponseFile);
end

if(exist('OFT_In_PreLinearisationCurve','var')==0)
    disp('using no linearization');
    OFT_In_PreLinearisationCurve='';
end

if(exist('OFT_In_IlluminationSpectrum','var')==0 || strcmp(OFT_In_IlluminationSpectrum,''))
    disp('using default daylight illumination'); 
    disp('D50');
    OFT_IlluminantSpectrum='D50'
else
    disp('using given illumination spectrum'); 
    disp(OFT_In_IlluminationSpectrum);
    OFT_IlluminantSpectrum=OFT_In_IlluminationSpectrum;
end

if(exist('OFT_In_SceneIllumination','var')==0 || strcmp(OFT_In_SceneIllumination,''))
    disp('using default scene daylight illumination'); 
    disp('D55');
    OFT_SceneIllumination='3400'
else
    disp('using given scene illumination spectrum'); 
    disp(OFT_In_SceneIllumination);
    OFT_SceneIllumination = OFT_In_SceneIllumination;
end

if(exist('OFT_In_PatchMeasurementFile','var')==0 || strcmp(OFT_In_PatchMeasurementFile,''))
    disp('using default patch set color checker'); 
    OFT_PatchMeasurementFile = 'Gretag Macbeth Color Checker'
else
    disp('using given patch set'); 
    disp(OFT_In_PatchMeasurementFile);
    OFT_PatchMeasurementFile = OFT_In_PatchMeasurementFile;
end

if(exist('OFT_In_NeutralsCompensation','var')==0)
    disp('using default neutral compensation');
    OFT_NeutralsCompensation=HDM_OFT_ColorNeutralCompensations.NoneType();  
else
    OFT_NeutralsCompensation=OFT_In_NeutralsCompensation;  
end

if(exist('OFT_In_ReferenceDomain','var')==0)
    disp('using default reference domain');   
    OFT_ReferenceDomain=HDM_OFT_IDT_ReferenceCamera.CIEType();   
else
    OFT_ReferenceDomain=OFT_In_ReferenceDomain;  
end

if(exist('OFT_In_StandardObserver','var')==0)
    disp('using default observer');   
    OFT_StandardObserver=HDM_OFT_CIEStandard.StandardObserver1931_2Degrees();   
else
    OFT_StandardObserver=OFT_In_StandardObserver;  
end

%% read and view linear input data

OFT_ImageOriginalOrg = HDM_OFT_ImageExportImport.ImportImage(OFT_In_CameraMeasurement, OFT_In_PreLinearisationCurve);
%OFT_ImageOriginalOrg = imresize(OFT_ImageOriginalOrg, 0.3);

if(isa(OFT_ImageOriginalOrg,'uint8'))
    OFT_ImageOriginalOrg=double(OFT_ImageOriginalOrg).*(1/(2^8-1));    
elseif(isa(OFT_ImageOriginalOrg,'uint16'))
   OFT_ImageOriginalOrg=double(OFT_ImageOriginalOrg).*(1/(2^16-1));  
   OFT_ImageOriginalOrg = imresize(OFT_ImageOriginalOrg, 0.4);
else
    OFT_ImageOriginalOrg = imresize(OFT_ImageOriginalOrg, 0.2);
end

[OFT_linUnprocessedInputImage_PatchLocations, OFT_linUnprocessedInputImage_PatchColours] = CCFind(OFT_ImageOriginalOrg);

if usejava('Desktop')
    visualizecc(OFT_ImageOriginalOrg, OFT_linUnprocessedInputImage_PatchLocations);
end

OFT_cameraImageOfTestChart=double(OFT_ImageOriginalOrg);

%% sRGB Curves:
sRGBDeLinearize = @(x)((x>0.0031308).*(1.055.*x.^(1/2.4)-0.055)+(x<=0.0031308).*12.92.*x);
sRGBLinearize   = @(x)((x>0.04045).*((x+0.055)./1.055).^2.4+(x<=0.04045).*(x./12.92));

figure('Name','original image - left/lin right/sRGB delinearized')
if usejava('Desktop')
	subplot(1,2,1),imshow(OFT_ImageOriginalOrg);
    subplot(1,2,2),imshow(sRGBDeLinearize(double(OFT_ImageOriginalOrg)));
end

% validate lin input
[OFT_linUnprocessedInputImage_PatchLocations, OFT_linUnprocessedInputImage_PatchColours] = CCFind(double(OFT_cameraImageOfTestChart));

if(isempty(OFT_linUnprocessedInputImage_PatchLocations))
    
    HDM_OFT_Utils.OFT_DispSubTitle('no patches found');
    
else
    
    % validate linearization
    l_GrayAxisData = OFT_linUnprocessedInputImage_PatchColours(:, 19:24);
    
    yR = polyfit(6:-1:1, l_GrayAxisData(1, :), 1);
    yG = polyfit(6:-1:1, l_GrayAxisData(2, :), 1);
    yB = polyfit(6:-1:1, l_GrayAxisData(3, :), 1);
    
    figure('Name','lin input gray axis patches values');
    if usejava('Desktop')
        plot(6:-1:1, l_GrayAxisData(1, :), 'r', ...
            6:-1:1, l_GrayAxisData(2, :), 'g', ...
            6:-1:1, l_GrayAxisData(3, :), 'b', ...
            6:-1:1, l_GrayAxisData(1, :), 'r+', ...
            6:-1:1, l_GrayAxisData(2, :), 'gx', ...
            6:-1:1, l_GrayAxisData(3, :), 'bo', ...
            6:-1:1, polyval(yR, 6:-1:1),'r--', ...
            6:-1:1, polyval(yG, 6:-1:1),'g--', ...
            6:-1:1, polyval(yB, 6:-1:1),'b--')
        xlabel('sample')
        ylabel('code value')
        legend({'r','g','b'})
        title(strcat('lin input gray axis patches values'));
    end

    % validate response
    
    % cam rgb
    l_cnames = {'dark skin', ...
                'light skin', ...
                'blue sky', ...
                'foliage', ...
                'blue flower', ...
                'bluish green', ...
                'orange', ...
                'purplish blue', ...
                'moderate red', ...
                'purple', ...
                'yellow green', ...
                'orange yellow', ...
                'blue', ...
                'green', ...
                'red', ...
                'yellow', ...
                'magenta', ...
                'cyan', ...
                'white 9.5 (.05 D)', ...
                'neutral 8 (.23 D)', ...
                'neutral 6.5 (.44 D)', ...
                'neutral 5 (.70 D)', ...
                'neutral 3.5 (1.05 D)', ...
                'black 2 (1.5 D)'};
      
    l_rnames = {'Red CAM','Green CAM','Blue CAM'};

    l_f = figure('Name','camera rgb values for patches');
    l_t = uitable(l_f ,'Data', OFT_linUnprocessedInputImage_PatchColours',...
                'ColumnName',l_rnames, ... 
                'Position',[0 0 100 100],...
                'RowName',l_cnames);
	l_t.Position(3) = l_t.Extent(3);
    l_t.Position(4) = l_t.Extent(4);
    
    OFT_PatchSetCameraTristimuliGreenNormalizedRatioFromCAM = OFT_linUnprocessedInputImage_PatchColours;
    for cur = 1 : size(OFT_linUnprocessedInputImage_PatchColours, 2)
        
        OFT_PatchSetCameraTristimuliGreenNormalizedRatioFromCAM(1, cur) = OFT_linUnprocessedInputImage_PatchColours(1, cur) / OFT_linUnprocessedInputImage_PatchColours(2, cur);
        OFT_PatchSetCameraTristimuliGreenNormalizedRatioFromCAM(2, cur) = OFT_linUnprocessedInputImage_PatchColours(2, cur) / OFT_linUnprocessedInputImage_PatchColours(2, cur);
        OFT_PatchSetCameraTristimuliGreenNormalizedRatioFromCAM(3, cur) = OFT_linUnprocessedInputImage_PatchColours(3, cur) / OFT_linUnprocessedInputImage_PatchColours(2, cur);
        
    end  
    
    l_f = figure('Name','camera rgb values green normalized ratio for patches by using lin input file');
    l_t = uitable(l_f ,'Data', OFT_PatchSetCameraTristimuliGreenNormalizedRatioFromCAM',...
                'ColumnName',l_rnames, ... 
                'Position',[0 0 100 100],...
                'RowName',l_cnames);
	l_t.Position(3) = l_t.Extent(3);
    l_t.Position(4) = l_t.Extent(4);
    
    % val from illumination/patch reflectances/camera response
    if(strfind(OFT_SpectralResponseFile, '.csv'))%if spectral based
        
        OFT_CameraSpectralResponse_1nm_CIE31Range = HDM_OFT_GetSpectralResponse(OFT_SpectralResponseFile);

        OFT_UsedIlluminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetIlluminantSpectrum(OFT_SceneIllumination);

        %%patches spectrum aquisition
        %10 nm resolution
        HDM_OFT_Utils.OFT_DispSubTitle('read patch spectra');
        OFT_PatchSet_SpectralCurve=HDM_OFT_PatchSet.GetPatchSpectra(OFT_PatchMeasurementFile);

        [OFT_PatchSetCameraTristimuli,OFT_PatchSetCameraTristimuli_ColorValueParts]=...
                    HDM_OFT_TristimuliCreator.CreateFromSpectrum(...
                            OFT_CameraSpectralResponse_1nm_CIE31Range,...
                            OFT_UsedIlluminant_Spectrum_1nm_CIE31Range,...
                            OFT_PatchSet_SpectralCurve);

        l_f = figure('Name','camera rgb values for patches by using illumination/reflectance/response');
        l_t = uitable(l_f ,'Data', OFT_PatchSetCameraTristimuli',...
                    'ColumnName',l_rnames, ... 
                    'Position',[0 0 100 100],...
                    'RowName',l_cnames);
        l_t.Position(3) = l_t.Extent(3);
        l_t.Position(4) = l_t.Extent(4);

        OFT_PatchSetCameraTristimuliGreenNormalizedRatioFromSpec = OFT_PatchSetCameraTristimuli;
        for cur = 1 : size(OFT_PatchSetCameraTristimuli, 2)

            OFT_PatchSetCameraTristimuliGreenNormalizedRatioFromSpec(1, cur) = OFT_PatchSetCameraTristimuli(1, cur) / OFT_PatchSetCameraTristimuli(2, cur);
            OFT_PatchSetCameraTristimuliGreenNormalizedRatioFromSpec(2, cur) = OFT_PatchSetCameraTristimuli(2, cur) / OFT_PatchSetCameraTristimuli(2, cur);
            OFT_PatchSetCameraTristimuliGreenNormalizedRatioFromSpec(3, cur) = OFT_PatchSetCameraTristimuli(3, cur) / OFT_PatchSetCameraTristimuli(2, cur);

        end  

        l_f = figure('Name','camera rgb values green normalized ratio for patches by using illumination/reflectance/response');
        l_t = uitable(l_f ,'Data', OFT_PatchSetCameraTristimuliGreenNormalizedRatioFromSpec',...
                    'ColumnName',l_rnames, ... 
                    'Position',[0 0 100 100],...
                    'RowName',l_cnames);
        l_t.Position(3) = l_t.Extent(3);
        l_t.Position(4) = l_t.Extent(4);


        l_ratioSpec2CAM = OFT_PatchSetCameraTristimuliGreenNormalizedRatioFromSpec ./ OFT_PatchSetCameraTristimuliGreenNormalizedRatioFromCAM;

        l_cnamesWithMeanAndSigma = l_cnames;
        l_cnamesWithMeanAndSigma{size(l_cnamesWithMeanAndSigma, 2) + 1} = 'Mean';
        l_cnamesWithMeanAndSigma{size(l_cnamesWithMeanAndSigma, 2) + 1} = 'Sigma';
        
        l_ratioSpec2CAM2View = l_ratioSpec2CAM';
        
        l_ratioSpec2CAM2View = ...
            [l_ratioSpec2CAM2View; ...
            mean(l_ratioSpec2CAM(1, :)), mean(l_ratioSpec2CAM(2, :)), mean(l_ratioSpec2CAM(3, :))];

        l_ratioSpec2CAM2View = ...
            [l_ratioSpec2CAM2View; ...
            std(l_ratioSpec2CAM(1, :)), std(l_ratioSpec2CAM(2, :)), std(l_ratioSpec2CAM(3, :))];
        
        l_f = figure('Name','spectral to CAM data ratio');
        l_t = uitable(l_f ,'Data', l_ratioSpec2CAM2View,...
                    'ColumnName',l_rnames, ... 
                    'Position',[0 0 100 100],...
                    'RowName',l_cnamesWithMeanAndSigma);
        l_t.Position(3) = l_t.Extent(3);
        l_t.Position(4) = l_t.Extent(4);
    
    end
    
end

%% process
% chapter 4.7.10 of AMPAS IDT standard: application of profile

%% gamma curve: non linear to linear

% for now nothing, cause input must be linear

%% target white point scale

OFT_TransformedImage2ViewbScaled = OFT_cameraImageOfTestChart;
OFT_TransformedImage2ViewbScaled(:,:,1) = (OFT_b(1)/min(OFT_b)) * OFT_cameraImageOfTestChart(:,:,1);
OFT_TransformedImage2ViewbScaled(:,:,2) = (OFT_b(2)/min(OFT_b)) * OFT_cameraImageOfTestChart(:,:,2);
OFT_TransformedImage2ViewbScaled(:,:,3) = (OFT_b(3)/min(OFT_b)) * OFT_cameraImageOfTestChart(:,:,3);

%% clip to max 1

OFT_TransformedImage2ViewbScaled(OFT_TransformedImage2ViewbScaled > 1.0) = 1.0;

%% apply matrix B, i.e. going from camera to ACES domain

OFT_TransformedImage2ViewBConv = OFT_TransformedImage2ViewbScaled;
OFT_TransformedImage2ViewBConv(:,:,1) = OFT_B(1,1)*OFT_TransformedImage2ViewbScaled(:,:,1) + OFT_B(1,2)*OFT_TransformedImage2ViewbScaled(:,:,2) + OFT_B(1,3)*OFT_TransformedImage2ViewbScaled(:,:,3);
OFT_TransformedImage2ViewBConv(:,:,2) = OFT_B(2,1)*OFT_TransformedImage2ViewbScaled(:,:,1) + OFT_B(2,2)*OFT_TransformedImage2ViewbScaled(:,:,2) + OFT_B(2,3)*OFT_TransformedImage2ViewbScaled(:,:,3);
OFT_TransformedImage2ViewBConv(:,:,3) = OFT_B(3,1)*OFT_TransformedImage2ViewbScaled(:,:,1) + OFT_B(3,2)*OFT_TransformedImage2ViewbScaled(:,:,2) + OFT_B(3,3)*OFT_TransformedImage2ViewbScaled(:,:,3);

figure('Name','b Scaled and  B Converted Image')
if usejava('Desktop')
	subplot(1,2,1),imshow(OFT_TransformedImage2ViewbScaled);
	subplot(1,2,2),imshow(OFT_TransformedImage2ViewBConv);
end

%% scale gray to 0.18

% only for test chart images possible
[OFT_cameraImageOfTestChart_PatchLocations,OFT_cameraImageOfTestChart_PatchColours] = CCFind(double(OFT_TransformedImage2ViewBConv));

if(isempty(OFT_cameraImageOfTestChart_PatchLocations))
    
    HDM_OFT_Utils.OFT_DispSubTitle('no patches found - no gray scaling done');
    
else
    
    l_18GrayVals = OFT_cameraImageOfTestChart_PatchColours(:, 22);
    
    l_facTo18Gray = 0.18 ./ l_18GrayVals;
    
    OFT_TransformedImage2ViewBConv(:, :, 1) = l_facTo18Gray(1) .* OFT_TransformedImage2ViewBConv(:, :, 1);
    OFT_TransformedImage2ViewBConv(:, :, 2) = l_facTo18Gray(2) .* OFT_TransformedImage2ViewBConv(:, :, 2);
    OFT_TransformedImage2ViewBConv(:, :, 3) = l_facTo18Gray(3) .* OFT_TransformedImage2ViewBConv(:, :, 3);
    
    if usejava('Desktop')
        visualizecc(OFT_TransformedImage2ViewBConv, OFT_cameraImageOfTestChart_PatchLocations);
    end    
       
    figure('Name','0.18 gray scaled, matrix transformed white normalized image - left/lin right/sRGB delinearized')
    if usejava('Desktop')
        subplot(1,2,1), imshow(OFT_TransformedImage2ViewBConv);
        subplot(1,2,2), imshow(sRGBDeLinearize(OFT_TransformedImage2ViewBConv));
    end
    
    [OFT_cameraImageOfTestChart_PatchLocations,OFT_cameraImageOfTestChart_PatchColours] = CCFind(double(OFT_TransformedImage2ViewBConv));
    
    if not(isempty(OFT_cameraImageOfTestChart_PatchLocations))
        %some bug to investigate: sometimes CCFind inverts the patch order
        if(OFT_cameraImageOfTestChart_PatchLocations(1,1) < OFT_cameraImageOfTestChart_PatchLocations(24,1))

            l_GrayAxisDataInACES = OFT_cameraImageOfTestChart_PatchColours(:, 19:24);

        else

            l_GrayAxisDataInACES = fliplr(OFT_cameraImageOfTestChart_PatchColours(:, 1:6));

        end

        yR = polyfit(6:-1:1, l_GrayAxisDataInACES(1, :), 1);
        yG = polyfit(6:-1:1, l_GrayAxisDataInACES(2, :), 1);
        yB = polyfit(6:-1:1, l_GrayAxisDataInACES(3, :), 1);

        figure('Name','gray axis patches values in ACES domain')
        if usejava('Desktop')
            plot(6:-1:1, l_GrayAxisDataInACES(1, :), 'r', ...
                6:-1:1, l_GrayAxisDataInACES(2, :), 'g', ...
                6:-1:1, l_GrayAxisDataInACES(3, :), 'b', ...
                6:-1:1, l_GrayAxisDataInACES(1, :), 'r+', ...
                6:-1:1, l_GrayAxisDataInACES(2, :), 'gx', ...
                6:-1:1, l_GrayAxisDataInACES(3, :), 'bo', ...
                6:-1:1, polyval(yR, 6:-1:1),'r--', ...
                6:-1:1, polyval(yG, 6:-1:1),'g--', ...
                6:-1:1, polyval(yB, 6:-1:1),'b--')
            xlabel('sample')
            ylabel('code value')
            legend({'r','g','b'})
            title(strcat('gray axis patches values in ACES domain'));
        end
    end
    
end

% now we might have a device independent color representation in aces color domain

%% convert from ACES to XYZ and to sRGB65 for viewing

[OFT_M,OFT_w]=HDM_OFT_IDT_ReferenceCamera.GetDefinition(OFT_ReferenceDomain);

OFT_Reference2StandardPrimaries=OFT_M;

OFT_MXYZ2sRGBD65=[3.2404542, -1.5371385, -0.4985314;
-0.9692660,  1.8760108,  0.0415560;
 0.0556434, -0.2040259,  1.0572252];

OFT_MOverall=OFT_MXYZ2sRGBD65 * OFT_Reference2StandardPrimaries;

OFT_TransformedImage2ViewBConv2XYZ2OutsRGB65 = OFT_TransformedImage2ViewBConv;
OFT_TransformedImage2ViewBConv2XYZ2OutsRGB65(:,:,1) = OFT_MOverall(1,1)*OFT_TransformedImage2ViewBConv(:,:,1) + OFT_MOverall(1,2)*OFT_TransformedImage2ViewBConv(:,:,2) + OFT_MOverall(1,3)*OFT_TransformedImage2ViewBConv(:,:,3);
OFT_TransformedImage2ViewBConv2XYZ2OutsRGB65(:,:,2) = OFT_MOverall(2,1)*OFT_TransformedImage2ViewBConv(:,:,1) + OFT_MOverall(2,2)*OFT_TransformedImage2ViewBConv(:,:,2) + OFT_MOverall(2,3)*OFT_TransformedImage2ViewBConv(:,:,3);
OFT_TransformedImage2ViewBConv2XYZ2OutsRGB65(:,:,3) = OFT_MOverall(3,1)*OFT_TransformedImage2ViewBConv(:,:,1) + OFT_MOverall(3,2)*OFT_TransformedImage2ViewBConv(:,:,2) + OFT_MOverall(3,3)*OFT_TransformedImage2ViewBConv(:,:,3);

figure('Name','matrix transformed white normalized image sRGB in D65 domain reillimunated - left/lin right/sRGB delinearized')
if usejava('Desktop')
	subplot(1,2,1),imshow(OFT_TransformedImage2ViewBConv2XYZ2OutsRGB65);
	subplot(1,2,2),imshow(sRGBDeLinearize(OFT_TransformedImage2ViewBConv2XYZ2OutsRGB65));
end

% compare xy and Lab steps

%% processed patches aces to xyz for patches

OFT_TransformedImage2View=OFT_TransformedImage2ViewBConv;

if(isempty(OFT_cameraImageOfTestChart_PatchLocations))
    HDM_OFT_Utils.OFT_DispSubTitle('no patches found no Lab evaluation');
    return ;
end

if usejava('Desktop')
    
    visualizecc(sRGBDeLinearize(OFT_TransformedImage2ViewBConv2XYZ2OutsRGB65) ,OFT_cameraImageOfTestChart_PatchLocations);
    
end

OFT_cameraImageOfTestChart_PatchColours_XYZConverted=OFT_Reference2StandardPrimaries * OFT_cameraImageOfTestChart_PatchColours;

%% processed patches XYZ and inverse chromatic adaptation

OFT_ChartImageTristimuli_NeutralsCompensated = OFT_cameraImageOfTestChart_PatchColours_XYZConverted;

% now we use the inverse for getting Illumination White from D60 based ACES
%OFT_ChartImageTristimuli_NeutralsCompensated=HDM_OFT_ColorNeutralCompensations.OFT_CompensateTristimuliForDifferentWhite...
%    (OFT_NeutralsCompensation, OFT_cameraImageOfTestChart_PatchColours_XYZConverted, OFT_w, OFT_Ww);

NumberOfPatches = size(OFT_ChartImageTristimuli_NeutralsCompensated,2);
OFT_ColorCheckerTristimuli_ColorValueParts = zeros(3,NumberOfPatches);

for patchIndex = 1:(NumberOfPatches) 
    cur = OFT_ChartImageTristimuli_NeutralsCompensated(:,patchIndex);
    l_sum = cur(1) + cur(2) + cur(3);
    OFT_ColorCheckerTristimuli_ColorValueParts(1,patchIndex) = cur(1) / l_sum;
    OFT_ColorCheckerTristimuli_ColorValueParts(2,patchIndex) = cur(2) / l_sum;
    OFT_ColorCheckerTristimuli_ColorValueParts(3,patchIndex) = cur(3) / l_sum;
end;

OFT_ChartImageTristimuli_NeutralsCompensatedNorm=OFT_ChartImageTristimuli_NeutralsCompensated;

OFT_ToLabWhite=OFT_w;

%% reference values from chart reflectances with target illumination and std observer

HDM_OFT_Utils.OFT_DispSubTitle('setup CIE standard observers curves');
OFT_CIEStandardObserver_SpectralCurves=HDM_OFT_CIEStandard.GetStandardObserverCurves(OFT_StandardObserver);

HDM_OFT_Utils.OFT_DispSubTitle('read patch spectra');
OFT_PatchSet_SpectralCurve=HDM_OFT_PatchSet.GetPatchSpectra(HDM_OFT_PatchSet.GretagMacbethColorChecker());

HDM_OFT_Utils.OFT_DispSubTitle('illuminat spectrum aquisition');
OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetIlluminantSpectrum(OFT_IlluminantSpectrum);

[OFT_ColorCheckerTristimuliReference,OFT_ColorCheckerTristimuli_ColorValueParts_Reference]=...
                HDM_OFT_TristimuliCreator.CreateFromSpectrum...
                    (OFT_CIEStandardObserver_SpectralCurves,...
                    OFT_Illuminant_Spectrum_1nm_CIE31Range,...
                    OFT_PatchSet_SpectralCurve);
               
%OFT_ColorCheckerTristimuliReference_NeutralsCompensated = OFT_ColorCheckerTristimuliReference;

OFT_Xw=trapz(OFT_CIEStandardObserver_SpectralCurves(2,:) .* OFT_Illuminant_Spectrum_1nm_CIE31Range(2,:));
OFT_Yw=trapz(OFT_CIEStandardObserver_SpectralCurves(3,:) .* OFT_Illuminant_Spectrum_1nm_CIE31Range(2,:));
OFT_Zw=trapz(OFT_CIEStandardObserver_SpectralCurves(4,:) .* OFT_Illuminant_Spectrum_1nm_CIE31Range(2,:));
OFT_WwUnscaled=[OFT_Xw,OFT_Yw,OFT_Zw]';
OFT_Ww=100*(OFT_WwUnscaled./OFT_WwUnscaled(2));

% now we use the inverse for getting Illumination White from D60 based ACES
OFT_ColorCheckerTristimuliReference_NeutralsCompensated=HDM_OFT_ColorNeutralCompensations.OFT_CompensateTristimuliForDifferentWhite...
    (OFT_NeutralsCompensation, OFT_ColorCheckerTristimuliReference, OFT_Ww, OFT_w);

NumberOfPatches = size(OFT_ColorCheckerTristimuliReference_NeutralsCompensated,2);
OFT_ColorCheckerTristimuli_ColorValueParts_Reference = zeros(3,NumberOfPatches);

for patchIndex = 1:(NumberOfPatches) 
    cur = OFT_ColorCheckerTristimuliReference_NeutralsCompensated(:,patchIndex);
    l_sum = cur(1) + cur(2) + cur(3);
    OFT_ColorCheckerTristimuli_ColorValueParts_Reference(1,patchIndex) = cur(1) / l_sum;
    OFT_ColorCheckerTristimuli_ColorValueParts_Reference(2,patchIndex) = cur(2) / l_sum;
    OFT_ColorCheckerTristimuli_ColorValueParts_Reference(3,patchIndex) = cur(3) / l_sum;
end;

OFT_ColorCheckerTristimuli_NeutralsCompensatedNorm = OFT_ColorCheckerTristimuliReference_NeutralsCompensated;

%% compute Lab for reference and image patches

OFT_ColorCheckerLab4Patches=HDM_OFT_ColorConversions.OFT_CIELab(OFT_ColorCheckerTristimuli_NeutralsCompensatedNorm, OFT_w);

% before we go to Lab domain for image XYZ: scale 0.18 gray patch to same Y as reference

l_XYZref4MidGray = OFT_ColorCheckerTristimuli_NeutralsCompensatedNorm(:, 22);

%some bug to investigate: sometimes CCFind inverts the patch order
if(OFT_cameraImageOfTestChart_PatchLocations(1,1) > OFT_cameraImageOfTestChart_PatchLocations(24,1))

    OFT_ChartImageTristimuli_NeutralsCompensatedNorm = fliplr(OFT_ChartImageTristimuli_NeutralsCompensatedNorm);

end

l_XYZimg4White = OFT_ChartImageTristimuli_NeutralsCompensatedNorm(:, 22);

l_scale4Ref = l_XYZref4MidGray ./ l_XYZimg4White;

OFT_ChartImageTristimuli_NeutralsCompensatedNorm(1, :) = OFT_ChartImageTristimuli_NeutralsCompensatedNorm(1, :).* l_scale4Ref(1,1);
OFT_ChartImageTristimuli_NeutralsCompensatedNorm(2, :) = OFT_ChartImageTristimuli_NeutralsCompensatedNorm(2, :).* l_scale4Ref(2,1);
OFT_ChartImageTristimuli_NeutralsCompensatedNorm(3, :) = OFT_ChartImageTristimuli_NeutralsCompensatedNorm(3, :).* l_scale4Ref(3,1);

OFT_ChartImageLab4Patches=HDM_OFT_ColorConversions.OFT_CIELab(OFT_ChartImageTristimuli_NeutralsCompensatedNorm, OFT_w);
    
%% compute delta E 

l_F = OFT_ColorCheckerLab4Patches - OFT_ChartImageLab4Patches;
l_deltaE4Patches = sqrt(sum(l_F.^2, 1));

l_deltaE2000 = deltaE2000(OFT_ColorCheckerLab4Patches', OFT_ChartImageLab4Patches');

l_rnames = {'L ref','a ref','b ref','L','a','b','delta E 2000'};

l_f = figure('Name','reference Lab values for patches (patch reflectances from BabelColor) and camera Lab values for patches for profiled image');
l_t = uitable(l_f ,'Data', [OFT_ColorCheckerLab4Patches', OFT_ChartImageLab4Patches', l_deltaE2000'],...
            'ColumnName',l_rnames, ... 
            'Position',[0 0 100 100],...
            'RowName',l_cnames);
l_t.Position(3) = l_t.Extent(3);
l_t.Position(4) = l_t.Extent(4);

OFT_MeanDeltaE4Patches=mean(l_deltaE2000);

deltaEStr=cellstr( num2str(l_deltaE2000','%4.3e\n') );

labels2 = cellstr( num2str([1:24]') );
labelsAnnot=labels2;
for cur=1:size(labelsAnnot)
    labelsAnnot(cur)=strcat(labelsAnnot(cur),'-- ');
end

oft_cld=strcat(labelsAnnot,deltaEStr);
oft_cld{end+1}='----------';
oft_cld{end+1}='mean delta E 2000:';
oft_cld{end+1}=num2str(OFT_MeanDeltaE4Patches);

%% getting spectral curve

OFT_CIE31_colorValueParts=csvread(strcat(OFT_Env.OFT_ConstraintsPath,'/cccie31.csv'));
OFT_CIE31_colorValueParts_x=[OFT_CIE31_colorValueParts(:,2)];
OFT_CIE31_colorValueParts_y=[OFT_CIE31_colorValueParts(:,3)];

if(strfind(OFT_SpectralResponseFile, '.csv'))%if spectral based
    
    %% compute saturated rgb xyz corner points by saturation

    l_perfectReflector = [OFT_PatchSet_SpectralCurve(1, :); ones(1, size(OFT_PatchSet_SpectralCurve, 2))];
    
    [l_Tristimuli4Corner_UsedIlluminant,l_Chromaticities4Corner_UsedIlluminant]=...
            HDM_OFT_TristimuliCreator.CreateFromSpectrum(...
                    OFT_CIEStandardObserver_SpectralCurves,...
                    OFT_UsedIlluminant_Spectrum_1nm_CIE31Range,...
                    OFT_CameraSpectralResponse_1nm_CIE31Range);        
                
    [l_Tristimuli4Corner_TargetIlluminant,l_Chromaticities4Corner_TargetIlluminant]=...
        HDM_OFT_TristimuliCreator.CreateFromSpectrum(...
                OFT_CIEStandardObserver_SpectralCurves,...
                OFT_Illuminant_Spectrum_1nm_CIE31Range,...
                OFT_CameraSpectralResponse_1nm_CIE31Range);      

    %% compute saturated rgb xyz corner points by matrix
    
    OFT_IDT_b;
    
    l_XYZ4R = OFT_Reference2StandardPrimaries * OFT_IDT_B * ([1; 0; 0] .* OFT_IDT_b);%' * l_M;
    l_xyz4R = l_XYZ4R / sum(l_XYZ4R(:));
    
    l_XYZ4G = OFT_Reference2StandardPrimaries * OFT_IDT_B * ([0; 1; 0] .* OFT_IDT_b);%' * l_M;
    l_xyz4G = l_XYZ4G / sum(l_XYZ4G(:));
    
    l_XYZ4B = OFT_Reference2StandardPrimaries * OFT_IDT_B * ([0; 0; 1] .* OFT_IDT_b);%' * l_M;
    l_xyz4B = l_XYZ4B / sum(l_XYZ4B(:));
    
    l_Chromaticities4Corner_MatrixBased = [l_xyz4R, l_xyz4G, l_xyz4B];

end

%% visualize

labelsRef = cellstr( num2str([1:24]') );

if(OFT_cameraImageOfTestChart_PatchLocations(1,1) < OFT_cameraImageOfTestChart_PatchLocations(24,1))
        
    labelsImg = labelsRef;

else

    labelsImg = cellstr( num2str([24:-1:1]') );

end

figure
plot([OFT_CIE31_colorValueParts_x;OFT_CIE31_colorValueParts_x(1)],[OFT_CIE31_colorValueParts_y;OFT_CIE31_colorValueParts_y(1)],'-',...
    OFT_ColorCheckerTristimuli_ColorValueParts(1,:), OFT_ColorCheckerTristimuli_ColorValueParts(2,:),'r+',...
    OFT_ColorCheckerTristimuli_ColorValueParts_Reference(1,:), OFT_ColorCheckerTristimuli_ColorValueParts_Reference(2,:),'bx', ...
    [l_Chromaticities4Corner_UsedIlluminant(1, :)'; l_Chromaticities4Corner_UsedIlluminant(1, 1)], ...
    [l_Chromaticities4Corner_UsedIlluminant(2, :)'; l_Chromaticities4Corner_UsedIlluminant(2, 1)],'r-', ...
    [l_Chromaticities4Corner_TargetIlluminant(1, :)'; l_Chromaticities4Corner_TargetIlluminant(1, 1)], ...
    [l_Chromaticities4Corner_TargetIlluminant(2, :)'; l_Chromaticities4Corner_TargetIlluminant(2, 1)],'g-', ...
    [l_Chromaticities4Corner_MatrixBased(1, :)'; l_Chromaticities4Corner_MatrixBased(1, 1)], ...
    [l_Chromaticities4Corner_MatrixBased(2, :)'; l_Chromaticities4Corner_MatrixBased(2, 1)],'b-')

text(OFT_ColorCheckerTristimuli_ColorValueParts_Reference(1,:), ...
    OFT_ColorCheckerTristimuli_ColorValueParts_Reference(2,:), labelsRef, 'VerticalAlignment','top', ...
    'HorizontalAlignment','left')

text(OFT_ColorCheckerTristimuli_ColorValueParts(1,:), ...
    OFT_ColorCheckerTristimuli_ColorValueParts(2,:), labelsImg, 'VerticalAlignment','bottom', ...
    'HorizontalAlignment','right')

xlabel('x')
ylabel('y')
legend({'','converted from estimated matrix transformed values','reference from babelcolor'})

annotation('textbox', [.15 .8, .1, .1],...
           'String', oft_cld,'FitHeightToText','on');

title('CIE31 x y color value parts');

OFT_TransformedImage2View = sRGBDeLinearize(OFT_TransformedImage2ViewBConv2XYZ2OutsRGB65);

if(exist('OFT_In_ReferenceImage2View','var'))
    % HDM_OFT_CompareIDTProfiledImage(OFT_In_ReferenceImage2View, OFT_TransformedImage2View);
end

HDM_OFT_Utils.OFT_DispTitle('evaluate transformed chart image succesfully finished');

end
