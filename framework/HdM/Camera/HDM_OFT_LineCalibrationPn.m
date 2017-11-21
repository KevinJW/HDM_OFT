function [out, o_lineR, o_lineG, o_lineB]=HDM_OFT_LineCalibrationPn...
    (i_SpectrometerMeasurement, i_CameraMeasurement, i_PreLinearisationCurve, i_maskImage)

OFT_Env=HDM_OFT_InitEnvironment();

HDM_OFT_Utils.OFT_DispTitle('start line calibration');

if(exist('i_SpectrometerMeasurement','var')==0)
    disp('using reference patch mesurements');
    i_SpectrometerMeasurement=strcat(OFT_IDTProcessorPath,'/cameraLineCalibrationReference/LineCali15042014.xls');
else
    disp(i_SpectrometerMeasurement);
end

if(exist('i_CameraMeasurement','var')==0)
    disp('using reference camera mesurements');
    i_CameraMeasurement=strcat(OFT_IDTProcessorPath,'/cameraLineCalibrationReference/Take021_Img0000050.TIF');
end

if(exist('i_PreLinearisationCurve','var')==0)
    disp('using no linearization');
    i_PreLinearisationCurve='';
end

%% read spectrometer data
HDM_OFT_Utils.OFT_DispSubTitle('read spectrometer data');

l_IntenityAgainstWavelength = HDM_OFT_SpectrumExportImport.ImportSpectrum(i_SpectrometerMeasurement);


%% find peaks
[peak_valueMat, peak_locationMat] = findpeaks(l_IntenityAgainstWavelength(2,:),...
    'minpeakheight',0.5*max(l_IntenityAgainstWavelength(2,:)),...
    'SORTSTR','descend');

[peak_value, peak_location] = HDM_OFT_findpeaks(l_IntenityAgainstWavelength(2,:),...
    9,1,100,...
    0.2 * max(l_IntenityAgainstWavelength(2,:)),...
    false);

[l_peaklocs, l_sorter] = sort(peak_location);
l_amps = zeros(size(peak_value));
l_npeaks = length(l_peaklocs);
for cur = 1:l_npeaks,
    l_amps(:,cur) = peak_value(:, l_sorter(cur));
end;

peak_value = l_amps;
peak_location = l_peaklocs;

%% average to close peaks

% [peak_location, peak_value] = HDM_OFT_AverageToClosePeaks(peak_location, peak_value, 2);

l_break = 0;

while l_break == 0
    
    [peak_location, peak_value] = ...
        HDM_OFT_AverageToClosePeaks(peak_location, peak_value, 6);

    l_ar = unique(uint16(peak_location));
       
    if size(l_ar, 2) < size(peak_location, 2)
        
        l_break = 0;
        
    else
        
        l_break = 1;
        
    end

end

%% create wavelength LUT
    
peak_location=int16(peak_location);

OFT_ReferencePeaksWaveLengths=zeros(1,size(peak_value,2));
OFT_ReferencePeaksWaveLengthsFiltered = [];
OFT_ReferencePeaksWaveLengthsValuesFiltered = [];
for peakIndex=1:size(peak_value,2)
    
    if l_IntenityAgainstWavelength(1,peak_location(peakIndex))==min(l_IntenityAgainstWavelength(1,:))
        continue;
    end
    
    if l_IntenityAgainstWavelength(1,peak_location(peakIndex))==max(l_IntenityAgainstWavelength(1,:))
        break;
    end
    
    peakEnv=[l_IntenityAgainstWavelength(1,peak_location(peakIndex)-1),l_IntenityAgainstWavelength(1,peak_location(peakIndex)),l_IntenityAgainstWavelength(1,peak_location(peakIndex)+1);...
        l_IntenityAgainstWavelength(2,peak_location(peakIndex)-1),peak_value(peakIndex),l_IntenityAgainstWavelength(2,peak_location(peakIndex)+1)];
    
    a=polyfit(peakEnv(1,:),peakEnv(2,:),2);    
    k = polyder(a);    
    r=roots(k);
    
    OFT_ReferencePeaksWaveLengths(1,peakIndex)=r;    
    polyval(a,r);
    
    if(r >425 && r < 650)
        
        OFT_ReferencePeaksWaveLengthsFiltered = [OFT_ReferencePeaksWaveLengthsFiltered, r];
        OFT_ReferencePeaksWaveLengthsValuesFiltered = [OFT_ReferencePeaksWaveLengthsValuesFiltered, peak_value(1, peakIndex)];
    
    end
    
end

OFT_ReferencePeaksWaveLengths = OFT_ReferencePeaksWaveLengthsFiltered;

figure('Name','Line Calibration Input')
subplot(2,2,1)
plot(l_IntenityAgainstWavelength(1,:),l_IntenityAgainstWavelength(2,:))

for cur = 1 : size(OFT_ReferencePeaksWaveLengths, 2)
    
    hold on
    
    plot([OFT_ReferencePeaksWaveLengths(1, cur) OFT_ReferencePeaksWaveLengths(1, cur)],ylim)
    
end

xlabel('wavelength in nm')
ylabel('intensity in W/(nm * m^2)')
title('line spectrum measured by spectrometer');


%% read camera line spectrum image
HDM_OFT_Utils.OFT_DispSubTitle('read camera line spectrum image');

%not aequidistant

%% find rect mask for spectrum region
l_maskImage = medfilt2(rgb2gray(HDM_OFT_ImageExportImport.ImportImage(i_maskImage)));
l_threshold = 0.005;%2 * (double(median(l_maskImage(20,:)))/double((1.0*2^16)));
l_maskImage1 = im2bw(l_maskImage, l_threshold);%//!!!

if usejava('Desktop')
	subplot(2,2,2),imshow(l_maskImage1);
end

%Find Plate6
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

%% read sopectrum image and masking

OFT_SpectrumImageFromFile=HDM_OFT_ImageExportImport.ImportImage(i_CameraMeasurement, i_PreLinearisationCurve);
[OFT_SpectrumImageFromFileHeight, OFT_SpectrumImageFromFileWidth,OFT_SpectrumImageFromFilefChannels]=size(OFT_SpectrumImageFromFile);

OFT_SpectrumImageFromFile(1:round(l_boundingBox4MaxArea(2)), :, :) = 0;
OFT_SpectrumImageFromFile(round(l_boundingBox4MaxArea(2) + l_boundingBox4MaxArea(4)) : size(OFT_SpectrumImageFromFile, 1), :, :) = 0;

OFT_SpectrumImageFromFile(:, 1:round(l_boundingBox4MaxArea(1)), :) = 0;
OFT_SpectrumImageFromFile(:, round(l_boundingBox4MaxArea(1) + l_boundingBox4MaxArea(3)) : size(OFT_SpectrumImageFromFile, 2), :) = 0;

l_subImage = imcrop(OFT_SpectrumImageFromFile, l_boundingBox4MaxArea);

if usejava('Desktop')

	subplot(2,2,3),imshow(OFT_SpectrumImageFromFile);
    
    subplot(2,2,2),imshow(l_subImage);

end

OFT_VCenterPos = round(l_boundingBox4MaxArea(1,2) + (l_boundingBox4MaxArea(1, 4)/2));

[RMax,RMaxInd]=max(OFT_SpectrumImageFromFile(OFT_VCenterPos,:,1));
[GMax,GMaxInd]=max(OFT_SpectrumImageFromFile(OFT_VCenterPos,:,2));

%if required rotate so that blue is always left
if(RMaxInd(1)<GMaxInd(1))
    OFT_SpectrumImage=flipdim(OFT_SpectrumImageFromFile,2);
else
    OFT_SpectrumImage=OFT_SpectrumImageFromFile;
end

OFT_SpectrumImageD = im2double(OFT_SpectrumImage);

OFT_SpectrumImageR = OFT_SpectrumImageD(:,:,1);
OFT_SpectrumImageG = OFT_SpectrumImageD(:,:,2);
OFT_SpectrumImageB = OFT_SpectrumImageD(:,:,3);

OFT_SpectrumImageR = OFT_SpectrumImageR / max(OFT_SpectrumImageR(:));
OFT_SpectrumImageG = OFT_SpectrumImageG / max(OFT_SpectrumImageG(:));
OFT_SpectrumImageB = OFT_SpectrumImageB / max(OFT_SpectrumImageB(:));

OFT_SpectrumImageGray = medfilt2(rgb2gray(OFT_SpectrumImage));

OFT_SpectrumImageGray1 = medfilt2((OFT_SpectrumImageR + OFT_SpectrumImageG + OFT_SpectrumImageB) / 3);

if usejava('Desktop')

	subplot(2,2,2),imshow(OFT_SpectrumImageGray1);
    
    subplot(2,2,3),imshow(OFT_SpectrumImage);
    
    subplot(2,2,4),imshow(OFT_SpectrumImageGray);

end

OFT_SpectrumROI_Top = round(l_boundingBox4MaxArea(2)) + 1;
OFT_SpectrumROI_Bottom = round(l_boundingBox4MaxArea(2) + l_boundingBox4MaxArea(4)) - 1;

OFT_Spectrum_MeanOfRows = zeros(1, size(OFT_SpectrumImageGray, 2));
OFT_Spectrum_MeanOfRowsR = zeros(1, size(OFT_SpectrumImageGray, 2));
OFT_Spectrum_MeanOfRowsG = zeros(1, size(OFT_SpectrumImageGray, 2));
OFT_Spectrum_MeanOfRowsB = zeros(1, size(OFT_SpectrumImageGray, 2));
OFT_Spectrum_TopRow = double(OFT_SpectrumImageGray(OFT_SpectrumROI_Top, :));
OFT_Spectrum_BottomRow = double(OFT_SpectrumImageGray(OFT_SpectrumROI_Bottom, :));
    
% norm
s1 = (OFT_Spectrum_TopRow - mean(OFT_Spectrum_TopRow)) / std(OFT_Spectrum_TopRow);
s2 = (OFT_Spectrum_BottomRow - mean(OFT_Spectrum_BottomRow)) / std(OFT_Spectrum_BottomRow);

% phase shift for rotation correction
%//!!! c = xcorr(s1, s2);                       %// Cross correlation
% lag = mod(find(c == max(c)), length(s2)); %// Find the position of the peak
% 
% %OFT_LineCaliTilt=(OFT_SpectrumImageWidth-lag)/(100)*360/(2*pi)/(OFT_SpectrumROI_Top-OFT_SpectrumROI_Bottom)
% OFT_LineCaliTilt=(lag/double((OFT_SpectrumROI_Bottom-OFT_SpectrumROI_Top)))*360.0/(2.0*pi);
% OFT_LineCaliTilt=0;%//!!!
% OFT_SpectrumImageGray=imrotate(OFT_SpectrumImageGray,OFT_LineCaliTilt,'crop');

% subplot(2,2,2),imshow(OFT_SpectrumImageIntensity);
%//!!!subplot(2,2,3),imshow(imcomplement(imrotate(OFT_SpectrumImageIntensity,OFT_LineCaliTilt)));

l_maxPerRow = zeros(5, OFT_SpectrumROI_Bottom - OFT_SpectrumROI_Top + 1);

l_maxPerRow(1,:) = 1 : size(l_maxPerRow, 2);

for R = OFT_SpectrumROI_Top : OFT_SpectrumROI_Bottom
    
    l_index = R - OFT_SpectrumROI_Top + 1;
    
    l_maxPerRow(1, l_index) = R;
    
    if R < size(OFT_SpectrumImageGray, 1)/2 - 50 || R > size(OFT_SpectrumImageGray, 1)/2 + 50
        continue;
    end
    
    l_maxPerRow(2, l_index) = max(OFT_SpectrumImageGray(R,:));
    
    cur = OFT_SpectrumImageGray(R,:);
    add = OFT_Spectrum_MeanOfRows(1,:);
    OFT_Spectrum_MeanOfRows = (double(add) + double(cur));
    
    l_maxPerRow(3, l_index) = max(OFT_SpectrumImageR(R,:));
    
    cur = OFT_SpectrumImageR(R,:);
    add = OFT_Spectrum_MeanOfRowsR(1,:);
    OFT_Spectrum_MeanOfRowsR = (double(add) + double(cur));
    
    l_maxPerRow(4, l_index) = max(OFT_SpectrumImageG(R,:));
    
    cur = OFT_SpectrumImageG(R,:);
    add = OFT_Spectrum_MeanOfRowsG(1,:);
    OFT_Spectrum_MeanOfRowsG = (double(add) + double(cur));
    
    l_maxPerRow(5, l_index) = max(OFT_SpectrumImageB(R,:));
    
    cur = OFT_SpectrumImageB(R,:);
    add = OFT_Spectrum_MeanOfRowsB(1,:);
    OFT_Spectrum_MeanOfRowsB = (double(add) + double(cur));
    
end

figure('Name','Max Distribution')
plot(l_maxPerRow(1,:),l_maxPerRow(2:end,:))

o_lineR = OFT_Spectrum_MeanOfRowsR;
o_lineG = OFT_Spectrum_MeanOfRowsG;
o_lineB = OFT_Spectrum_MeanOfRowsB;

OFT_Spectrum_MeanOfRows=OFT_Spectrum_MeanOfRows/cast((OFT_SpectrumROI_Bottom-OFT_SpectrumROI_Top),'like',OFT_Spectrum_MeanOfRows);
OFT_Spectrum_MeanOfRows=OFT_Spectrum_MeanOfRows/max(OFT_Spectrum_MeanOfRows(1,:));
OFT_Spectrum_MeanOfRowsFlipped=fliplr(OFT_Spectrum_MeanOfRows);

subplot(2,2,4)
plot(s1,'g');
hold on
plot(s2,'b');
hold on
plot(OFT_Spectrum_MeanOfRows(1,:)/max(OFT_Spectrum_MeanOfRows(1,:)),'r');
legend({'Top Row:(i-mean(I))/\sigma (I)','Bottom Row:(i-mean(I))/\sigma (I)','All Rows: normalized mean value'})
xlabel('horicontal pixel index')
ylabel('pixel value')
title('detected lines (vertical gray bars)');

threshold=0.1;%//!!! 0.2 for BL75 0.12 for zup 40
% matlab based find peaks from certain package
%[OFT_SpectrumImage_peak_valueMat, OFT_SpectrumImage_peak_locationMat] = ...
%    findpeaks(OFT_Spectrum_MeanOfRows,'MINPEAKHEIGHT',threshold * max(OFT_Spectrum_MeanOfRows),...
%    'SORTSTR','descend');%...
%    %,'THRESHOLD',0.01*(min(OFT_Spectrum_MeanOfRows)+1))
    
% [OFT_SpectrumImage_peak_valueFlippedMat, OFT_SpectrumImage_peak_locationFlippedMat] = ...
%     findpeaks(OFT_Spectrum_MeanOfRowsFlipped,'MINPEAKHEIGHT',threshold * max(OFT_Spectrum_MeanOfRowsFlipped),...
%     'SORTSTR','descend');%...
%     %,'THRESHOLD',0.01*(min(OFT_Spectrum_MeanOfRows)+1))
    
% own find peaks    
[OFT_SpectrumImage_peak_value, OFT_SpectrumImage_peak_location] = HDM_OFT_findpeaks(OFT_Spectrum_MeanOfRows,...
    2 * size(OFT_ReferencePeaksWaveLengths, 2),10,100,...
    threshold * max(OFT_Spectrum_MeanOfRows),...
    false);

[l_peaklocs, l_sorter] = sort(OFT_SpectrumImage_peak_location);
l_amps = zeros(size(OFT_SpectrumImage_peak_value));
l_npeaks = length(l_peaklocs);
for cur = 1:l_npeaks,
    l_amps(:,cur) = OFT_SpectrumImage_peak_value(:, l_sorter(cur));
end;

OFT_SpectrumImage_peak_value = l_amps;
OFT_SpectrumImage_peak_location = l_peaklocs;

l_break = 0;

l_max  = l_peaklocs(1, size(l_peaklocs, 2));
l_min = l_peaklocs(1, 1);
l_delta = l_max - l_min;

l_averageWidth = l_delta / (2.0 * size(l_peaklocs, 2)); % unit pixel

while l_break == 0
    
    [OFT_SpectrumImage_peak_location, OFT_SpectrumImage_peak_value] = ...
        HDM_OFT_AverageToClosePeaks(OFT_SpectrumImage_peak_location, OFT_SpectrumImage_peak_value, l_averageWidth);

    l_ar = unique(uint16(OFT_SpectrumImage_peak_location));
       
    if size(l_ar, 2) < size(OFT_SpectrumImage_peak_location, 2)
        
        l_break = 0;
        
    else
        
        l_break = 1;
        
    end

end

[OFT_SpectrumImage_peak_valueFlipped, OFT_SpectrumImage_peak_locationFlipped] = HDM_OFT_findpeaks(OFT_Spectrum_MeanOfRowsFlipped,...
    2 * size(OFT_ReferencePeaksWaveLengths, 2),10,100,...
    threshold * max(OFT_Spectrum_MeanOfRowsFlipped),...
    false);

[l_peaklocs, l_sorter] = sort(OFT_SpectrumImage_peak_locationFlipped);
l_amps = zeros(size(OFT_SpectrumImage_peak_valueFlipped));
l_npeaks = length(l_peaklocs);
for cur = 1:l_npeaks,
    l_amps(:,cur) = OFT_SpectrumImage_peak_valueFlipped(:, l_sorter(cur));
end;

OFT_SpectrumImage_peak_valueFlipped = l_amps;
OFT_SpectrumImage_peak_locationFlipped = l_peaklocs;

l_break = 0;

while l_break == 0
    
    [OFT_SpectrumImage_peak_locationFlipped, OFT_SpectrumImage_peak_valueFlipped] = ...
        HDM_OFT_AverageToClosePeaks(OFT_SpectrumImage_peak_locationFlipped, OFT_SpectrumImage_peak_valueFlipped, l_averageWidth);

    l_ar = unique(uint16(OFT_SpectrumImage_peak_locationFlipped));
       
    if size(l_ar, 2) < size(OFT_SpectrumImage_peak_locationFlipped, 2)
        
        l_break = 0;
        
    else
        
        l_break = 1;
        
    end

end

l_OFT_Spectrum_MeanOfRowsFiltered = sgolayfilt(OFT_Spectrum_MeanOfRows, 3, 11);
    
width=size(OFT_SpectrumImageGray,2);    
OFT_SpectrumImage_peak_locationFlipped=(-1*OFT_SpectrumImage_peak_locationFlipped)+width;

OFT_SpectrumImage_peak_locationFlippedReverse = fliplr(OFT_SpectrumImage_peak_locationFlipped);

l_maxIndex = size(OFT_SpectrumImage_peak_location,2);

if l_maxIndex > size(OFT_SpectrumImage_peak_locationFlippedReverse,2);
    
    l_maxIndex = size(OFT_SpectrumImage_peak_locationFlippedReverse,2);
    
end

OFT_SpectrumImage_peak_value = OFT_SpectrumImage_peak_value(1,1:l_maxIndex);

OFT_SpectrumImage_peak_location = 0.5 * (OFT_SpectrumImage_peak_location(1,1:l_maxIndex) + OFT_SpectrumImage_peak_locationFlippedReverse(1,1:l_maxIndex));
l_break = 0;

while l_break == 0
    
    [OFT_SpectrumImage_peak_location, OFT_SpectrumImage_peak_value] = ...
        HDM_OFT_AverageToClosePeaks(OFT_SpectrumImage_peak_location, OFT_SpectrumImage_peak_value, l_averageWidth);

    l_ar = unique(uint16(OFT_SpectrumImage_peak_location));
       
    if size(l_ar, 2) < size(OFT_SpectrumImage_peak_location, 2)
        
        l_break = 0;
        
    else
        
        l_break = 1;
        
    end

end

x = sort(OFT_SpectrumImage_peak_location(1,:));
l_intensity4x = zeros(1, size(OFT_SpectrumImage_peak_value, 2));

for cur = 1 : size(x, 2)
    
    pos = x(1, cur);
    l_lookup = find(OFT_SpectrumImage_peak_location == pos);
    l_intensity4x(cur) = OFT_SpectrumImage_peak_value(1, l_lookup(1, 1));
      
end

[l_amps, l_sorter] = sort(OFT_SpectrumImage_peak_value);
l_peaklocs = zeros(size(OFT_SpectrumImage_peak_location));
l_npeaks = length(l_amps);
for cur = 1:l_npeaks,
    l_peaklocs(:,cur) = OFT_SpectrumImage_peak_location(:, l_sorter(cur));
end;

OFT_SpectrumImage_peak_value_IntensitySorted = fliplr(l_amps);
OFT_SpectrumImage_peak_location_IntensitySorted = fliplr(l_peaklocs);

if OFT_SpectrumImage_peak_location_IntensitySorted(1,1) > OFT_SpectrumImage_peak_location_IntensitySorted(1, 2)
    
    l_lGreenPeekPos = find(l_intensity4x == OFT_SpectrumImage_peak_value_IntensitySorted(1, 2));
    l_rRedPeekPos = find(l_intensity4x == OFT_SpectrumImage_peak_value_IntensitySorted(1, 1));
    
else
    
    l_lGreenPeekPos = find(l_intensity4x == OFT_SpectrumImage_peak_value_IntensitySorted(1, 1));
    l_rRedPeekPos = find(l_intensity4x == OFT_SpectrumImage_peak_value_IntensitySorted(1, 2));    
    
end 

v = sort(OFT_ReferencePeaksWaveLengths(1,:));
l_intensity4v = zeros(1, size(OFT_ReferencePeaksWaveLengthsValuesFiltered, 2));

for cur = 1 : size(v, 2)
    
    pos = v(1, cur);
    l_lookup = find(OFT_ReferencePeaksWaveLengths == pos);
    l_intensity4v(cur) = OFT_ReferencePeaksWaveLengthsValuesFiltered(1, l_lookup);
      
end

[l_amps, l_sorter] = sort(OFT_ReferencePeaksWaveLengthsValuesFiltered);
l_peaklocs = zeros(size(OFT_ReferencePeaksWaveLengths));
l_npeaks = length(l_amps);
for cur = 1:l_npeaks,
    l_peaklocs(:,cur) = OFT_ReferencePeaksWaveLengths(:, l_sorter(cur));
end;

OFT_ReferencePeaksWaveLengthsValuesFiltered_IntensitySorted = fliplr(l_amps);
OFT_ReferencePeaksWaveLengths_IntensitySorted = fliplr(l_peaklocs);

if OFT_ReferencePeaksWaveLengths_IntensitySorted(1,1) > OFT_ReferencePeaksWaveLengths_IntensitySorted(1, 2)
    
    l_lGreenPeekPosRef = find(l_intensity4v == OFT_ReferencePeaksWaveLengthsValuesFiltered_IntensitySorted(1, 2));
    l_rRedPeekPosRef = find(l_intensity4v == OFT_ReferencePeaksWaveLengthsValuesFiltered_IntensitySorted(1, 1));
    
else
    
    l_lGreenPeekPosRef = find(l_intensity4v == OFT_ReferencePeaksWaveLengthsValuesFiltered_IntensitySorted(1, 1));
    l_rRedPeekPosRef = find(l_intensity4v == OFT_ReferencePeaksWaveLengthsValuesFiltered_IntensitySorted(1, 2));    
    
end 

x1 = zeros(1, min(l_lGreenPeekPos, l_lGreenPeekPosRef) + 1 + min(size(x, 2) - l_rRedPeekPos, size(v, 2) - l_rRedPeekPosRef));
v1 = x1;

for cur = 1 : min(l_lGreenPeekPos, l_lGreenPeekPosRef) - 1 
    
    x1(min(l_lGreenPeekPos, l_lGreenPeekPosRef) - cur) = x(l_lGreenPeekPos - cur);
    v1(min(l_lGreenPeekPos, l_lGreenPeekPosRef) - cur) = v(l_lGreenPeekPosRef - cur);
    
end

x1(min(l_lGreenPeekPos, l_lGreenPeekPosRef)) = x(l_lGreenPeekPos);
v1(min(l_lGreenPeekPos, l_lGreenPeekPosRef)) = v(l_lGreenPeekPosRef);

x1(min(l_lGreenPeekPos, l_lGreenPeekPosRef) + 1) = x(l_rRedPeekPos);
v1(min(l_lGreenPeekPos, l_lGreenPeekPosRef) + 1) = v(l_rRedPeekPosRef);

l_cnt = 1;
for cur = min(l_lGreenPeekPos, l_lGreenPeekPosRef) + 2 : size(x1, 2)
    
    x1(cur) = x(l_rRedPeekPos + l_cnt);
    v1(cur) = v(l_rRedPeekPosRef + l_cnt);
    l_cnt = l_cnt + 1;
    
end

x = x1;
v = v1;

% //!!! we do not take the righter outmost redline for now
% //!!! to be discussed by adding a red laser line
l_indexOfLastUsedRedLine = size(v, 2);
OFT_SpectrumImagePixelColumnIndex=interp1(x(1:l_indexOfLastUsedRedLine), v(1:l_indexOfLastUsedRedLine), 1 : 1 : size(OFT_SpectrumImageGray, 2),'linear','extrap');%//!!!'linear',

l_WaveLengthOfHorizontalPixelCoordinatePolynomial = polyfit(x(1:l_indexOfLastUsedRedLine),v(1:l_indexOfLastUsedRedLine), 2);

l_WaveLengthOfHorizontalPixelCoordinate = polyval(l_WaveLengthOfHorizontalPixelCoordinatePolynomial, 1 : size(OFT_SpectrumImageGray, 2));

% OFT_SpectrumImagePixelColumnIndex = l_WaveLengthOfHorizontalPixelCoordinate;

%%found lines in image

figure('Name','Lines Detected in Image')

plot(OFT_Spectrum_MeanOfRows)

hold on

plot(l_OFT_Spectrum_MeanOfRowsFiltered)

for cur = 1 : size(x, 2)
    
    hold on
    
    plot([x(1, cur) x(1, cur)],ylim)
        
end

xlabel('wavelength in nm')
ylabel('normalized intensity')
title('average line spectrum measured by camera');

%%result

figure('Name','Line Calibration Result')
plot(l_IntenityAgainstWavelength(1,:),l_IntenityAgainstWavelength(2,:)/max(l_IntenityAgainstWavelength(2,:)))
hold on
plot(OFT_SpectrumImagePixelColumnIndex,OFT_Spectrum_MeanOfRows(1,:)/max(OFT_Spectrum_MeanOfRows(1,:)),'r');

for cur = 1 : size(OFT_ReferencePeaksWaveLengths, 2)
    
    hold on
    plot([OFT_ReferencePeaksWaveLengths(1, cur) OFT_ReferencePeaksWaveLengths(1, cur)],ylim)
    
end

xlim([350 850])

legend({'normalized reference spectrometer data','normalized mean value of camera spectrum'})
xlabel('mapped hypothetical wavelength in nm')
ylabel('relative intensity')
title('found and mapped two dominant lines (vertical gray bars)');

out=OFT_SpectrumImagePixelColumnIndex;

HDM_OFT_Utils.OFT_DispTitle('line calibration succesfully finished');

end

%//!!!
%MEDFILT1       One-dimensional median filter
%
%       y = MEDFILT(x)
%       y = MEDFILT(x, w)
%
%       median filter the signal with window of width W (default is 5).

% Copyright (C) 1995-2009, by Peter I. Corke
%
% This file is part of The Machine Vision Toolbox for Matlab (MVTB).
% 
% MVTB is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MVTB is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with MVTB.  If not, see <http://www.gnu.org/licenses/>.
function m = medfilt1(s, w)
        if nargin == 1,
                w = 5;
        end
        
        s = s(:)';
        w2 = floor(w/2);
        w = 2*w2 + 1;

        n = length(s);
        m = zeros(w,n+w-1);
        s0 = s(1); sl = s(n);

        for i=0:(w-1),
                m(i+1,:) = [s0*ones(1,i) s sl*ones(1,w-i-1)];
        end
        m = median(m);
        m = m(w2+1:w2+n);
end