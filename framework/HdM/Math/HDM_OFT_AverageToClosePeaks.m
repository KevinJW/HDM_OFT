function [o_Peak_location, o_Peak_value] = HDM_OFT_AverageToClosePeaks(i_Peak_location, i_Peak_value, i_width)

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
