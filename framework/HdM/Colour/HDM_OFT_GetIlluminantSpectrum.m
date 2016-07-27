function out=HDM_OFT_GetIlluminantSpectrum(OFT_IlluminantSpectrum)

if(isempty(strfind(OFT_IlluminantSpectrum,'.')))
    
    val = str2num(OFT_IlluminantSpectrum);
    
    if(isempty(val))
    
        switch OFT_IlluminantSpectrum

            case 'A'
                OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetBlackBodyRadiatorIllumination(2856);

            case '2800'
                OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetBlackBodyRadiatorIllumination(2800);
            case '3000'
                OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetBlackBodyRadiatorIllumination(3000);
            case '3050'
                OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetBlackBodyRadiatorIllumination(3050);
            case 'ISOTungsten'
                OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetBlackBodyRadiatorIllumination(3050);
            case '3200'
                OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetBlackBodyRadiatorIllumination(3200);
            case '3400'
                OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetBlackBodyRadiatorIllumination(3400);

            case 'C'
                OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetStandardCIllumination();%//!!!transpose avoid

            case 'D50'
                OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetStandardDIllumination(5000);
            case 'D55'
                OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetStandardDIllumination(5500);
            case 'D60'
                OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetStandardDIllumination(6000);
            case 'D65'
                OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetStandardDIllumination(6504);
            case 'D75'
                OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetStandardDIllumination(7500);

            case 'E'

                l_wavelengths = 360:1:830;
                l_intensity = ones(1, size(l_wavelengths, 2));
                OFT_Illuminant_Spectrum_1nm_CIE31Range = [l_wavelengths; l_intensity];

            otherwise
        end
        
    else
        
        if val < 4500
            
            OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetBlackBodyRadiatorIllumination(val);
            
        else
            
            OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetStandardDIllumination(val);
            
        end
        
    end
else
    %OFT_Illuminant_Spectrum_1nm_CIE31Range = HDM_OFT_GetCIERange1nmSpectrumFromSpectralDataFile(OFT_IlluminantSpectrum);
    OFT_Illuminant_Spectrum_1nm_CIE31Range = HDM_OFT_SpectrumExportImport.ImportSpectrum(OFT_IlluminantSpectrum);
end

out=OFT_Illuminant_Spectrum_1nm_CIE31Range;

end

