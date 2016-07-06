function OFT_Out_Planck_Spectrum_1nm_CIE31Range = HDM_OFT_GetBlackBodyRadiatorIllumination(OFT_In_Temperature)

OFT_Env=HDM_OFT_InitEnvironment();

HDM_OFT_Utils.OFT_DispTitle('Planck Black Body Radiator');

if(exist('OFT_In_Temperature','var')==0)   
    disp('using default temperature');
    OFT_BlackBodyTemperature=3200
else    
    default = 0;    
    disp('using given temperature');
    OFT_BlackBodyTemperature=OFT_In_Temperature
end

%% black body spectrum
HDM_OFT_Utils.OFT_DispSubTitle('Black Body Spectrum');

l_c = 3*10^8;           % speed of light in vaccum
l_h = 6.625*10.^-34;    % Planck constant 
k = 1.38*10.^-23;       % Boltzmann constant

l_lambda=(0.36:0.001:0.830).*1e-6; 

l_intensity = (2 * l_h * l_c * l_c) ./ ((l_lambda.^5) .* (exp((l_h .* l_c) ./ (k .* OFT_BlackBodyTemperature .* l_lambda)) - 1));
l_intensity = l_intensity ./ max(l_intensity(:));

OFT_Out_Planck_Spectrum_1nm_CIE31Range=[360:1:830 ; l_intensity];

if(exist('default','var')==0)
    figure
    plot(OFT_Out_Planck_Spectrum_1nm_CIE31Range(1,:),OFT_Out_Planck_Spectrum_1nm_CIE31Range(2,:))
    xlabel('wavelength in nm')
    ylabel('relative power')
    title(strcat('Black Body Radiation for ',num2str(OFT_BlackBodyTemperature),' K'));
end

HDM_OFT_Utils.OFT_DispTitle('Planck Black Body Radiator finished');

end

