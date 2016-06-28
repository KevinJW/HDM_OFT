function out=HDM_OFT_LightCalibration...
    (OFT_In_Pixel2WavelengthLookUp, ...
    OFT_In_SpectrometerMeasurement, OFT_In_CameraMeasurement,...
    OFT_In_PreLinearisationCurve,...
    OFT_In_Sensor, OFT_In_FocalLength)

OFT_Env=HDM_OFT_InitEnvironment();

HDM_OFT_Utils.OFT_DispTitle('start light calibration');

if(exist('OFT_In_SpectrometerMeasurement','var')==0)
    disp('using reference patch mesurements');
    OFT_In_SpectrometerMeasurement=strcat(OFT_IDTProcessorPath,'/../../cameraLightCalibrationReference/LightCali15042014.xls');
else
    disp(OFT_In_SpectrometerMeasurement);
end

if(exist('OFT_In_CameraMeasurement','var')==0)
    OFT_In_CameraMeasurement=strcat(OFT_IDTProcessorPath,'/../../cameraLightCalibrationReference/Take022_Img0000050.TIF');
    disp('using reference camera mesurements');
else
    disp(OFT_In_CameraMeasurement);  
end

if(exist('OFT_In_PreLinearisationCurve','var')==0)
    disp('using no linearization');
    OFT_In_PreLinearisationCurve='';
end

%% read spectrometer data
HDM_OFT_Utils.OFT_DispSubTitle('read spectrometer data');

l_IntenityAgainstWavelength = HDM_OFT_SpectrumExportImport.ImportSpectrum(OFT_In_SpectrometerMeasurement);

% //!!! to be discussed for noisy uprtec data 
window = 10;
h = ones(window,1)/window;
l_IntenityAgainstWavelength2 = l_IntenityAgainstWavelength;

for cur = 5 : size(l_IntenityAgainstWavelength2, 2) - 5
    
    l_nonZeroCnt = 0;
    l_sum = 0;
    
    for innerCur = cur - 4 : cur + 4
        
        if(l_IntenityAgainstWavelength(2, innerCur) > 0)
            
            l_nonZeroCnt = l_nonZeroCnt + 1;
            
            l_sum = l_sum + l_IntenityAgainstWavelength(2, innerCur);
            
        end
        
    end
    
    if (l_sum>0 && l_nonZeroCnt>0)
        
        l_IntenityAgainstWavelength2(2, cur) = l_sum / l_nonZeroCnt;
    
    end
         
end

l_IntenityAgainstWavelength = l_IntenityAgainstWavelength2;

figure
subplot(2,2,1)
plot(l_IntenityAgainstWavelength(1,:),l_IntenityAgainstWavelength(2,:), ...
    l_IntenityAgainstWavelength2(1,:),l_IntenityAgainstWavelength2(2,:))
xlabel('wavelength in nm')
ylabel('intensity in W/(nm * m^2)')
title('tungsten light spectrum measured by reference spectrometer');

% //!!! l_IntenityAgainstWavelength = l_IntenityAgainstWavelength2;

%% read camera line spectrum image
HDM_OFT_Utils.OFT_DispSubTitle('read camera line spectrum image');

%not aequidistant

OFT_SpectrumImage = HDM_OFT_ImageExportImport.ImportImage(OFT_In_CameraMeasurement, OFT_In_PreLinearisationCurve);

%% find rect mask for spectrum region
l_maskImage = medfilt2(rgb2gray(OFT_SpectrumImage));
l_maskImage1 = im2bw(l_maskImage, 0.05);

%Find Plate
[lab, n] = bwlabel(l_maskImage1);

regions = regionprops(lab, 'All');
regionsCount = size(regions, 1) ;

l_area = 0;
l_boundingBox4MaxArea = [];

for i = 1:regionsCount
    
    region = regions(i);
    
    if(l_area < region.Area)
        
        l_boundingBox4MaxArea = region.BoundingBox;
        l_area = region.Area;
        
    end

end

%% change mask dimensions

l_boundingBox4MaxArea(1) = l_boundingBox4MaxArea(1) - 0.5 * l_boundingBox4MaxArea(3);
l_boundingBox4MaxArea(3) = 2 * l_boundingBox4MaxArea(3);

l_boundingBox4MaxArea(2) = l_boundingBox4MaxArea(2) + 0.1 * l_boundingBox4MaxArea(4);
l_boundingBox4MaxArea(4) = l_boundingBox4MaxArea(4) - 0.2 * l_boundingBox4MaxArea(4);

%% apply mask

OFT_SpectrumImage(1:round(l_boundingBox4MaxArea(2)), :, :) = 0;
OFT_SpectrumImage(round(l_boundingBox4MaxArea(2) + l_boundingBox4MaxArea(4)) : size(OFT_SpectrumImage, 1), :, :) = 0;

OFT_SpectrumImage(:, 1:round(l_boundingBox4MaxArea(1)), :) = 0;
OFT_SpectrumImage(:, round(l_boundingBox4MaxArea(1) + l_boundingBox4MaxArea(3)) : size(OFT_SpectrumImage, 2), :) = 0;

%% median filtering

OFT_SpectrumImage(:,:,1) = medfilt2(OFT_SpectrumImage(:,:,1));
OFT_SpectrumImage(:,:,2) = medfilt2(OFT_SpectrumImage(:,:,2));
OFT_SpectrumImage(:,:,3) = medfilt2(OFT_SpectrumImage(:,:,3));

%//!!! OFT_SpectrumImage=HDM_OFT_Cos4Correction(OFT_SpectrumImage,OFT_In_FocalLength,OFT_In_Sensor);%//!!!
[OFT_SpectrumImageHeight, OFT_SpectrumImageWidth,OFT_SpectrumImageNofChannels]=size(OFT_SpectrumImage);

OFT_VCenterPos = round(l_boundingBox4MaxArea(1,2) + (l_boundingBox4MaxArea(1, 4)/2));

[RMax,RMaxInd]=max(OFT_SpectrumImage(OFT_VCenterPos,:,1));
[GMax,GMaxInd]=max(OFT_SpectrumImage(OFT_VCenterPos,:,2));

%% if required rotate so that blue is always left

%if required rotate so that blue is always left
if(RMaxInd(1)<GMaxInd(1))
    OFT_SpectrumImage = flipdim(OFT_SpectrumImage,2);
    
    l_boundingBox4MaxArea(1) = size(OFT_SpectrumImage, 2) - round(l_boundingBox4MaxArea(1) + l_boundingBox4MaxArea(3));
    
end

if usejava('Desktop')
	subplot(2,2,2),imshow(OFT_SpectrumImage);
    
    l_subImage = imcrop(OFT_SpectrumImage, l_boundingBox4MaxArea);   
    subplot(2,2,3),imshow(l_subImage);
    
end


OFT_Center = OFT_VCenterPos;

OFT_SpectrumROI_Top = round(l_boundingBox4MaxArea(2));
OFT_SpectrumROI_Bottom = round(l_boundingBox4MaxArea(2) + l_boundingBox4MaxArea(4));

OFT_SpectrumImageGray=rgb2gray(OFT_SpectrumImage);
OFT_Spectrum_MeanOfRows=zeros(1,OFT_SpectrumImageWidth,OFT_SpectrumImageNofChannels);

%//!!!OFT_SpectrumImage=imrotate(OFT_SpectrumImage,1,'crop');

for R=OFT_SpectrumROI_Top:OFT_SpectrumROI_Bottom
    cur=OFT_SpectrumImage(R,:,:);
    add=OFT_Spectrum_MeanOfRows(1,:,:);
    OFT_Spectrum_MeanOfRows=(double(add)+double(cur));
end

OFT_Spectrum_MeanOfRows=OFT_Spectrum_MeanOfRows/cast((OFT_SpectrumROI_Bottom-OFT_SpectrumROI_Top),'like',OFT_Spectrum_MeanOfRows);

OFT_SpectrumImagePixelColumnIndex=1:OFT_SpectrumImageWidth;

subplot(2,2,4)
plot...
    (OFT_SpectrumImagePixelColumnIndex,OFT_Spectrum_MeanOfRows(1,:,3),...
    OFT_SpectrumImagePixelColumnIndex,OFT_Spectrum_MeanOfRows(1,:,2),...
    OFT_SpectrumImagePixelColumnIndex,OFT_Spectrum_MeanOfRows(1,:,1));
legend({'b','g','r'})
xlabel('horicontal pixel index')
ylabel('mean pixel value of all lines')

OFT_first= round(l_boundingBox4MaxArea(1, 1));
OFT_last = round(l_boundingBox4MaxArea(1, 1) + l_boundingBox4MaxArea(1, 3));

OFT_In_Pixel2WavelengthLookUp=OFT_In_Pixel2WavelengthLookUp(OFT_first:OFT_last);
OFT_Spectrum_MeanOfRows=OFT_Spectrum_MeanOfRows(1,OFT_first:OFT_last,:);

OFT_Spectrum_MeanOfRowsRed_nmSamples=...
    interp1(OFT_In_Pixel2WavelengthLookUp,OFT_Spectrum_MeanOfRows(1,:,1),...
    min(l_IntenityAgainstWavelength(1,:)):l_IntenityAgainstWavelength(1,2)-l_IntenityAgainstWavelength(1,1):max(l_IntenityAgainstWavelength(1,:)),...
    'linear');

OFT_Spectrum_MeanOfRowsGreen_nmSamples=...
    interp1(OFT_In_Pixel2WavelengthLookUp,OFT_Spectrum_MeanOfRows(1,:,2),...
    min(l_IntenityAgainstWavelength(1,:)):l_IntenityAgainstWavelength(1,2)-l_IntenityAgainstWavelength(1,1):max(l_IntenityAgainstWavelength(1,:)),...
    'linear');

OFT_Spectrum_MeanOfRowsBlue_nmSamples=...
    interp1(OFT_In_Pixel2WavelengthLookUp,OFT_Spectrum_MeanOfRows(1,:,3),...
    min(l_IntenityAgainstWavelength(1,:)):l_IntenityAgainstWavelength(1,2)-l_IntenityAgainstWavelength(1,1):max(l_IntenityAgainstWavelength(1,:)),...
    'linear');

OFT_Spectrum_MeanOfRowsRed_nmSamplesCam = OFT_Spectrum_MeanOfRowsRed_nmSamples;
OFT_Spectrum_MeanOfRowsGreen_nmSamplesCam = OFT_Spectrum_MeanOfRowsGreen_nmSamples;
OFT_Spectrum_MeanOfRowsBlue_nmSamplesCam = OFT_Spectrum_MeanOfRowsBlue_nmSamples;

l_maxGreen = max(OFT_Spectrum_MeanOfRowsGreen_nmSamplesCam);
OFT_Spectrum_MeanOfRowsRed_nmSamplesCam = OFT_Spectrum_MeanOfRowsRed_nmSamplesCam ./ l_maxGreen;
OFT_Spectrum_MeanOfRowsGreen_nmSamplesCam = OFT_Spectrum_MeanOfRowsGreen_nmSamplesCam ./ l_maxGreen;
OFT_Spectrum_MeanOfRowsBlue_nmSamplesCam = OFT_Spectrum_MeanOfRowsBlue_nmSamplesCam ./ l_maxGreen;

OFT_NormalizedLight = l_IntenityAgainstWavelength(2,:)/max(l_IntenityAgainstWavelength(2,:));
OFT_NormalizedLightCorrection = OFT_NormalizedLight;
for cur = 1 : size(OFT_NormalizedLightCorrection, 2)
    
    if (OFT_NormalizedLight(1, cur) > 0)
        
        OFT_NormalizedLightCorrection(1, cur) = 1/OFT_NormalizedLight(1, cur);
        
    end

end
OFT_NormalizedLightCorrection=OFT_NormalizedLightCorrection ./ max(OFT_NormalizedLightCorrection);

OFT_Spectrum_MeanOfRowsRed_nmSamples=OFT_Spectrum_MeanOfRowsRed_nmSamples.*OFT_NormalizedLightCorrection;
OFT_Spectrum_MeanOfRowsGreen_nmSamples=OFT_Spectrum_MeanOfRowsGreen_nmSamples.*OFT_NormalizedLightCorrection;
OFT_Spectrum_MeanOfRowsBlue_nmSamples=OFT_Spectrum_MeanOfRowsBlue_nmSamples.*OFT_NormalizedLightCorrection;

l_maxGreen = max(OFT_Spectrum_MeanOfRowsGreen_nmSamples);
OFT_Spectrum_MeanOfRowsRed_nmSamples = OFT_Spectrum_MeanOfRowsRed_nmSamples ./ l_maxGreen;
OFT_Spectrum_MeanOfRowsGreen_nmSamples = OFT_Spectrum_MeanOfRowsGreen_nmSamples ./ l_maxGreen;
OFT_Spectrum_MeanOfRowsBlue_nmSamples = OFT_Spectrum_MeanOfRowsBlue_nmSamples ./ l_maxGreen;

figure ('Name','light and light correction')
plot...
    (l_IntenityAgainstWavelength(1,:),OFT_NormalizedLight(1,:),...
    l_IntenityAgainstWavelength(1,:),OFT_NormalizedLightCorrection(1,:));
legend({'normalized SPD of Tungsten light', 'light correction curve'})
xlabel('wavelength in nm')
ylabel('relative radiance')
grid on
grid minor

figure ('Name','response uncorrected and light corrected')
plot...
    (l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsBlue_nmSamples,...
    l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsGreen_nmSamples,...
    l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsRed_nmSamples,...
    l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsBlue_nmSamplesCam,...
    l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsGreen_nmSamplesCam,...
    l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsRed_nmSamplesCam);
legend({'b light corrected','g light corrected','r light corrected','b from CAM','g from CAM','r from CAM'})
xlabel('wavelength in nm')
ylabel('pixel value of all lines')
grid on
grid minor

[oft_ndata, oft_text, oft_alldata] =xlsread(strcat(OFT_Env.OFT_ConstraintsPath,'/RuledGratingEfficiency.xlsx'));
OFT_GridTransmission=oft_ndata;

l_lowerValidWavelength = 360;
l_upperValidWavelength = 830;

for cur = 2 : size(l_IntenityAgainstWavelength, 2) - 1
    
    if (l_IntenityAgainstWavelength(2, cur) > 0 && l_IntenityAgainstWavelength(2, cur - 1) == 0)
        
        l_lowerValidWavelength = l_IntenityAgainstWavelength(1, cur);
        
    end
    
    if (l_IntenityAgainstWavelength(2, cur) == 0 && l_IntenityAgainstWavelength(2, cur - 1) > 0)
        
        l_upperValidWavelength = l_IntenityAgainstWavelength(1, cur);
        
    end
    
end

l_nmBase = 360 : 1 : 830;
OFT_GridTransmissionNormIntensity=interp1(OFT_GridTransmission(1,:), OFT_GridTransmission(2,:), 360 : 1 : 830, 'pchip', 'extrap');

l_norm = 1;
for cur = 1 : size(l_nmBase, 2)
    
    if(l_nmBase(1, cur) > l_upperValidWavelength)
        
        l_norm = OFT_GridTransmissionNormIntensity(1, cur-1);
        break;
        
    end
    
end    

OFT_GridTransmissionNormIntensity=1/l_norm*OFT_GridTransmissionNormIntensity;

OFT_GridTransmission_Correction=1./OFT_GridTransmissionNormIntensity;

l_norm = 1;
for cur = 1 : size(l_nmBase, 2)
    
    if(l_nmBase(1, cur) > l_upperValidWavelength)
        
        l_norm = OFT_GridTransmission_Correction(1, cur-1);
        break;
        
    end
    
end   

OFT_GridTransmission_Correction = OFT_GridTransmission_Correction ./ l_norm;

OFT_Spectrum_MeanOfRowsRed_nmSamplesLight=OFT_Spectrum_MeanOfRowsRed_nmSamples;
OFT_Spectrum_MeanOfRowsGreen_nmSamplesLight=OFT_Spectrum_MeanOfRowsGreen_nmSamples;
OFT_Spectrum_MeanOfRowsBlue_nmSamplesLight=OFT_Spectrum_MeanOfRowsBlue_nmSamples;

OFT_Spectrum_MeanOfRowsRed_nmSamples=OFT_GridTransmission_Correction.*OFT_Spectrum_MeanOfRowsRed_nmSamples;%//!!! 1.3
OFT_Spectrum_MeanOfRowsGreen_nmSamples=OFT_GridTransmission_Correction.*OFT_Spectrum_MeanOfRowsGreen_nmSamples;
OFT_Spectrum_MeanOfRowsBlue_nmSamples=OFT_GridTransmission_Correction.*OFT_Spectrum_MeanOfRowsBlue_nmSamples;%//!!! 0.9

l_maxGreen = max(OFT_Spectrum_MeanOfRowsGreen_nmSamples);
OFT_Spectrum_MeanOfRowsRed_nmSamples = OFT_Spectrum_MeanOfRowsRed_nmSamples ./ l_maxGreen;
OFT_Spectrum_MeanOfRowsGreen_nmSamples = OFT_Spectrum_MeanOfRowsGreen_nmSamples ./ l_maxGreen;
OFT_Spectrum_MeanOfRowsBlue_nmSamples = OFT_Spectrum_MeanOfRowsBlue_nmSamples ./ l_maxGreen;


figure ('Name','normalized grating efficiency and grating correction')
plot...
    (l_IntenityAgainstWavelength(1,:),OFT_GridTransmissionNormIntensity(1,:),...
    l_IntenityAgainstWavelength(1,:),OFT_GridTransmission_Correction(1,:))
legend({'normalized grating efficiency', 'grating correction'})
xlabel('wavelength in nm')
ylabel('normalized efficiency and correction')
grid on
grid minor

l_lightAndGratingCorrection = OFT_NormalizedLightCorrection .* OFT_GridTransmission_Correction;
l_lightAndGratingCorrection = l_lightAndGratingCorrection ./ max(l_lightAndGratingCorrection);

figure ('Name','grating and light correction')
plot...
    (l_IntenityAgainstWavelength(1,:), OFT_NormalizedLightCorrection, ...
    l_IntenityAgainstWavelength(1,:), OFT_GridTransmission_Correction, ...
    l_IntenityAgainstWavelength(1,:), l_lightAndGratingCorrection)
legend({'light correction', 'grating correction', 'light and grating correction'})
xlabel('wavelength in nm')
ylabel('correction factor')
grid on
grid minor

figure ('Name','response light corrected and grating corrected')
plot...
    (l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsBlue_nmSamples,...
    l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsGreen_nmSamples,...
    l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsRed_nmSamples,...
    l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsBlue_nmSamplesLight,...
    l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsGreen_nmSamplesLight,...
    l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsRed_nmSamplesLight)
legend({'b grating corrected','g grating corrected','r grating corrected','b light corrected','g light corrected','r light corrected'})
xlabel('wavelength in nm')
ylabel('green normalized pixel value of all lines (grating corrected)')
grid on
grid minor

OFT_Out_SpectralResponse=[l_IntenityAgainstWavelength(1,:);...
    OFT_Spectrum_MeanOfRowsRed_nmSamples/max(OFT_Spectrum_MeanOfRowsGreen_nmSamples);...
    OFT_Spectrum_MeanOfRowsGreen_nmSamples/max(OFT_Spectrum_MeanOfRowsGreen_nmSamples);...
    OFT_Spectrum_MeanOfRowsBlue_nmSamples/max(OFT_Spectrum_MeanOfRowsGreen_nmSamples)];

figure ('Name','final response')
plot...
    (l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsBlue_nmSamples/max(OFT_Spectrum_MeanOfRowsGreen_nmSamples),...
    l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsGreen_nmSamples/max(OFT_Spectrum_MeanOfRowsGreen_nmSamples),...
    l_IntenityAgainstWavelength(1,:),OFT_Spectrum_MeanOfRowsRed_nmSamples/max(OFT_Spectrum_MeanOfRowsGreen_nmSamples))
legend({'b','g','r','grid efficency'})
xlabel('wavelength in nm')
ylabel('green normalized pixel value of all lines (phase grid transmission corrected)')
grid on
grid minor

% for test: show one line for noise evaluation
% plot...
% 	(OFT_SpectrumImagePixelColumnIndex,OFT_SpectrumImage(OFT_SpectrumROI_Top+(OFT_SpectrumROI_Bottom-OFT_SpectrumROI_Top)/2,1:OFT_SpectrumImageWidth,1),'r',...
%     OFT_SpectrumImagePixelColumnIndex,OFT_SpectrumImage(OFT_SpectrumROI_Top+(OFT_SpectrumROI_Bottom-OFT_SpectrumROI_Top)/2,1:OFT_SpectrumImageWidth,2),'g',...
%     OFT_SpectrumImagePixelColumnIndex,OFT_SpectrumImage(OFT_SpectrumROI_Top+(OFT_SpectrumROI_Bottom-OFT_SpectrumROI_Top)/2,1:OFT_SpectrumImageWidth,3),'b');
% 
%     %(OFT_SpectrumImagePixelColumnIndex,OFT_SpectrumImage(OFT_SpectrumROI_Top+(OFT_SpectrumROI_Bottom-OFT_SpectrumROI_Top)/2,1:OFT_SpectrumImageWidth,1),OFT_SpectrumImagePixelColumnIndex,OFT_SpectrumImage(OFT_SpectrumROI_Top+(OFT_SpectrumROI_Bottom-OFT_SpectrumROI_Top)/2,1:OFT_SpectrumImageWidth,2),OFT_SpectrumImagePixelColumnIndex,OFT_SpectrumImage(OFT_SpectrumROI_Top+(OFT_SpectrumROI_Bottom-OFT_SpectrumROI_Top)/2,1:OFT_SpectrumImageWidth,3));
% 
% xlabel('horicontal pixel index')
% ylabel('red pixel value of center line')

out=OFT_Out_SpectralResponse;

HDM_OFT_Utils.OFT_DispTitle('light calibration succesfully finished');

end

%//!!! open source functions if avoiding matlab image toolbox usage

function [rout,g,b] = imresize(varargin)
%IMRESIZE Resize image.
%   B = IMRESIZE(A,M,'method') returns an image matrix that is 
%   M times larger (or smaller) than the image A.  The image B
%   is computed by interpolating using the method in the string
%   'method'.  Possible methods are 'nearest' (nearest neighbor),
%   'bilinear' (binlinear interpolation), or 'bicubic' (bicubic 
%   interpolation). B = IMRESIZE(A,M) uses 'nearest' as the 
%   default interpolation scheme.
%
%   B = IMRESIZE(A,[MROWS NCOLS],'method') returns a matrix of 
%   size MROWS-by-NCOLS.
%
%   RGB1 = IMRESIZE(RGB,...) resizes the RGB truecolor image 
%   stored in the 3-D array RGB, and returns a 3-D array (RGB1).
%
%   When the image size is being reduced, IMRESIZE lowpass filters
%   the image before interpolating to avoid aliasing. By default,
%   this filter is designed using FIR1, but can be specified 
%   using IMRESIZE(...,'method',H).  The default filter is 11-by-11.
%   IMRESIZE(...,'method',N) uses an N-by-N filter.
%   IMRESIZE(...,'method',0) turns off the filtering.
%   Unless a filter H is specified, IMRESIZE will not filter
%   when 'nearest' is used.
%   
%   See also IMZOOM, FIR1, INTERP2.

%   Grandfathered Syntaxes:
%
%   [R1,G1,B1] = IMRESIZE(R,G,B,M,'method') or 
%   [R1,G1,B1] = IMRESIZE(R,G,B,[MROWS NCOLS],'method') resizes
%   the RGB image in the matrices R,G,B.  'bilinear' is the
%   default interpolation method.

%   Clay M. Thompson 7-7-92
%   Copyright (c) 1992 by The MathWorks, Inc.
%   $Revision: 5.4 $  $Date: 1996/10/16 20:33:27 $

[A,m,method,classIn,h] = parse_inputs(varargin{:});

threeD = (ndims(A)==3); % Determine if input includes a 3-D array

if threeD,
   r = resizeImage(A(:,:,1),m,method,h);
   g = resizeImage(A(:,:,2),m,method,h);
   b = resizeImage(A(:,:,3),m,method,h);
   if nargout==0, 
      imshow(r,g,b);
      return;
   elseif nargout==1,
      if strcmp(classIn,'uint8');
         rout = repmat(uint8(0),[size(r),3]);
         rout(:,:,1) = uint8(round(r*255));
         rout(:,:,2) = uint8(round(g*255));
         rout(:,:,3) = uint8(round(b*255));
      else
         rout = zeros([size(r),3]);
         rout(:,:,1) = r;
         rout(:,:,2) = g;
         rout(:,:,3) = b;
      end
   else % nargout==3
      if strcmp(classIn,'uint8')
         rout = uint8(round(r*255)); 
         g = uint8(round(g*255)); 
         b = uint8(round(b*255)); 
      else
         rout = r;        % g,b are already defined correctly above
      end
   end
else 
   r = resizeImage(A,m,method,h);
   if nargout==0,
      imshow(r);
      return;
   end
   if strcmp(classIn,'uint8')
      r = uint8(round(r*255)); 
   end
   rout = r;
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function: resizeImage
%

function b = resizeImage(A,m,method,h)
% Inputs:
%         A       Input Image
%         m       resizing factor or 1-by-2 size vector
%         method  'nearest','bilinear', or 'bicubic'
%         h       the anti-aliasing filter to use.
%                 if h is zero, don't filter
%                 if h is an integer, design and use a filter of size h
%                 if h is empty, use default filter

if prod(size(m))==1,
   bsize = floor(m*size(A));
else
   bsize = m;
end

if any(size(bsize)~=[1 2]),
   error('M must be either a scalar multiplier or a 1-by-2 size vector.');
end

% values in bsize must be at least 1.
bsize = max(bsize, 1);

if (any((bsize < 4) & (bsize < size(A))) & ~strcmp(method, 'nea'))
   fprintf('Input is too small for bilinear or bicubic method;\n');
   fprintf('using nearest-neighbor method instead.\n');
   method = 'nea';
end

if isempty(h),
   nn = 11; % Default filter size
else
   if prod(size(h))==1, 
      nn = h; h = []; 
   else 
      nn = 0;
   end
end

[m,n] = size(A);

if nn>0 & method(1)=='b',  % Design anti-aliasing filter if necessary
   if bsize(1)1 | length(h2)>1, h = h1'*h2; else h = []; end
   if length(h1)>1 | length(h2)>1, 
      a = filter2(h1',filter2(h2,A)); 
   else 
      a = A; 
   end
elseif method(1)=='b' & (prod(size(h)) > 1),
   a = filter2(h,A);
else
   a = A;
end

if method(1)=='n', % Nearest neighbor interpolation
   dx = n/bsize(2); dy = m/bsize(1); 
   uu = (dx/2+.5):dx:n+.5; vv = (dy/2+.5):dy:m+.5;
elseif all(method == 'bil') | all(method == 'bic'),
   uu = 1:(n-1)/(bsize(2)-1):n; vv = 1:(m-1)/(bsize(1)-1):m;
else
   error(['Unknown interpolation method: ',method]);
end

%
% Interpolate in blocks
%
nu = length(uu); nv = length(vv);
blk = bestblk([nv nu]);
nblks = floor([nv nu]./blk); nrem = [nv nu] - nblks.*blk;
mblocks = nblks(1); nblocks = nblks(2);
mb = blk(1); nb = blk(2);

rows = 1:blk(1); b = zeros(nv,nu);
for i=0:mblocks,
   if i==mblocks, rows = (1:nrem(1)); end
   for j=0:nblocks,
      if j==0, cols = 1:blk(2); elseif j==nblocks, cols=(1:nrem(2)); end
      if ~isempty(rows) & ~isempty(cols)
         [u,v] = meshgrid(uu(j*nb+cols),vv(i*mb+rows));
         % Interpolate points
         if method(1) == 'n', % Nearest neighbor interpolation
            b(i*mb+rows,j*nb+cols) = interp2(a,u,v,'*nearest');
         elseif all(method == 'bil'), % Bilinear interpolation
            b(i*mb+rows,j*nb+cols) = interp2(a,u,v,'*linear');
         elseif all(method == 'bic'), % Bicubic interpolation
            b(i*mb+rows,j*nb+cols) = interp2(a,u,v,'*cubic');
         end
      end
   end
end

if nargout==0,
   if isgray(b), imshow(b,size(colormap,1)), else imshow(b), end
   return
end

if isgray(A)   % This should always be true
   b = max(0,min(b,1));  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function: parse_inputs
%

function [A,m,method,classIn,h] = parse_inputs(varargin)
% Outputs:  A       the input image
%           m       the resize scaling factor or the new size
%           method  interpolation method (nearest,bilinear,bicubic)
%           class   storage class of A
%           h       if 0, skip filtering; if non-zero scalar, use filter 
%                   of size h; otherwise h is the anti-aliasing filter.

switch nargin
case 2,                        % imresize(A,m)
   A = varargin{1};
   m = varargin{2};
   method = 'nearest';
   classIn = class(A);
   h = [];
case 3,                        % imresize(A,m,method)
   A = varargin{1};
   m = varargin{2};
   method = varargin{3};
   classIn = class(A);
   h = [];
case 4,
   if isstr(varargin{3})       % imresize(A,m,method,h)
      A = varargin{1};
      m = varargin{2};
      method = varargin{3};
      classIn = class(A);
      h = varargin{4};
   else                        % imresize(r,g,b,m)
      for i=1:3
         if isa(varargin{i},'uint8')
            error('Please use 3-d RGB array syntax with uint8 image data');
         end
      end
      A = zeros([size(varargin{1}),3]);
      A(:,:,1) = varargin{1};
      A(:,:,2) = varargin{2};
      A(:,:,3) = varargin{3};
      m = varargin{4};
      method = 'nearest';
      classIn = class(A);
      h = [];
   end
case 5,                        % imresize(r,g,b,m,'method')
   for i=1:3
      if isa(varargin{i},'uint8')
         error('Please use 3-d RGB array syntax with uint8 image data');
      end
   end
   A = zeros([size(varargin{1}),3]);
   A(:,:,1) = varargin{1};
   A(:,:,2) = varargin{2};
   A(:,:,3) = varargin{3};
   m = varargin{4};
   method = varargin{5};
   classIn = class(A);
   h = [];
case 6,                        % imresize(r,g,b,m,'method',h)
   for i=1:3
      if isa(varargin{i},'uint8')
         error('Please use 3-d RGB array syntax with uint8 image data');
      end
   end
   A = zeros([size(varargin{1}),3]);
   A(:,:,1) = varargin{1};
   A(:,:,2) = varargin{2};
   A(:,:,3) = varargin{3};
   m = varargin{4};
   method = varargin{5};
   classIn = class(A);
   h = varargin{6};
otherwise,
   error('Invalid input arguments.');
end

if isa(A, 'uint8'),     % Convert A to Double grayscale for filtering & interpolation
   A = double(A)/255;
end

method = [lower(method),'   ']; % Protect against short method
method = method(1:3);
end



function [mb,nb] = bestblk(siz,k)
%BESTBLK Best block size for block processing.
%	BLK = BESTBLK([M N],K) returns the 1-by-2 block size BLK
%	closest to but smaller than K-by-K for block processing.
%
%	[MB,NB] = BESTBLK([M N],K) returns the best block size
%	as the two scalars MB and NB.
%
%	[...] = BESTBLK([M N]) returns the best block size smaller
%	than 100-by-100.
%
%	BESTBLK returns the M or N when they are already smaller
%	than K.
%
%	See also BLKPROC, SIZE.

%	Clay M. Thompson
%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.6 $  $Date: 1994/03/04 19:54:04 $

if nargin==1, k = 100; end % Default block size

%
% Find possible factors of siz that make good blocks
%

% Define acceptable block sizes
m = floor(k):-1:floor(min(ceil(siz(1)/10),k/2));
n = floor(k):-1:floor(min(ceil(siz(2)/10),k/2));

% Choose that largest acceptable block that has the minimum padding.
[dum,ndx] = min(ceil(siz(1)./m).*m-siz(1)); blk(1) = m(ndx);
[dum,ndx] = min(ceil(siz(2)./n).*n-siz(2)); blk(2) = n(ndx);

if nargout==2,
  mb = blk(1); nb = blk(2);
else
  mb = blk;
end

end


function y = isgray(x)
%ISGRAY True for intensity images.
%	ISGRAY(A) returns 1 if A is an intensity image and 0 otherwise.
%	An intensity image contains values between 0.0 and 1.0.
%
%	See also ISIND, ISBW.

%	Clay M. Thompson 2-25-93
%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.4 $  $Date: 1993/08/18 03:11:32 $

y = min(min(x))>=0 & max(max(x))<=1;
end