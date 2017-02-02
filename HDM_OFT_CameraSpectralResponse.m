function o_cameraSpectralResponse = HDM_OFT_CameraSpectralResponse(OFT_In_IDTTaskData)

OFT_Env=HDM_OFT_InitEnvironment();

HDM_OFT_Utils.OFT_DispTitle('camera spectral response estimation');

%% pixel to wavelength mapping using line light

[l_Pixel2WavelengthLookUp_HigherOrderPolynomBased, l_LineR, l_LineG, l_LineB] = HDM_OFT_LineCalibrationPn...
    (OFT_In_IDTTaskData.SpectralResponse_In_LineCalibrationSpectrum, OFT_In_IDTTaskData.SpectralResponse_In_LineCalibrationImage,...
    OFT_In_IDTTaskData.PreLinearisation_Out_LinCurve, OFT_In_IDTTaskData.SpectralResponse_In_LightCalibrationImage);

%% normalized camera response estimation using tungsten light

l_CameraResponse = HDM_OFT_LightCalibrationPn...
    (l_Pixel2WavelengthLookUp_HigherOrderPolynomBased,...
    OFT_In_IDTTaskData.SpectralResponse_In_LightCalibrationSpectrum, OFT_In_IDTTaskData.SpectralResponse_In_LightCalibrationImage,...
    OFT_In_IDTTaskData.PreLinearisation_Out_LinCurve, OFT_In_IDTTaskData.Device_In_Sensor, OFT_In_IDTTaskData.Device_In_FocalLength);

%% error correction

o_cameraSpectralResponse = HDM_OFT_CameraResponseCorrection(l_CameraResponse, OFT_In_IDTTaskData.SpectralResponse_In_LineCalibrationSpectrum, ...
    l_LineR, l_LineG, l_LineB, ...
    l_Pixel2WavelengthLookUp_HigherOrderPolynomBased);

% o_cameraSpectralResponse = l_CameraResponse;

HDM_OFT_Utils.OFT_DispTitle('camera spectral response estimation succesfully finished');

end
