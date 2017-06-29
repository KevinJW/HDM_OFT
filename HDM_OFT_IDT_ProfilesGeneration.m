function OFT_IDTFiles=HDM_OFT_IDT_ProfilesGeneration...
    (IDTTaskData)

if isempty(IDTTaskData.SpectralResponse_In_LightCalibrationSpectrum)        %chart based
    OFT_Illuminations={IDTTaskData.IDTCreationConstraints_In_WhitePoint};
else                                                                        %spectral based therefore we use the line and light spectra too
    OFT_Illuminations={IDTTaskData.IDTCreationConstraints_In_WhitePoint};%, ...
                        %IDTTaskData.SpectralResponse_In_LineCalibrationSpectrum, IDTTaskData.SpectralResponse_In_LightCalibrationSpectrum};
end

%adding the AMPAS specified light sources
if strcmp(OFT_Illuminations(1), 'D55')
    
    OFT_Illuminations{size(OFT_Illuminations, 2) + 1} = '3050';
    
else
    
    OFT_Illuminations{size(OFT_Illuminations, 2) + 1} = 'D55'; 
    OFT_Illuminations{size(OFT_Illuminations, 2) + 1} = '3050';

end

OFT_IDTFiles=cell(size(OFT_Illuminations));

HDM_OFT_Utils.OFT_DispTitle('unify patch sets');

if not(isempty(IDTTaskData.IDTCreationConstraints_In_AdditionalPatchSets))
    
    l_objectReflectances = HDM_OFT_PatchSet.GetPatchSpectra(IDTTaskData.IDTCreationConstraints_In_PatchSet);
    
    for cur = 1 : size(IDTTaskData.IDTCreationConstraints_In_AdditionalPatchSets, 1)
        
        l_curObjectReflectances = HDM_OFT_PatchSet.GetPatchSpectra(IDTTaskData.IDTCreationConstraints_In_AdditionalPatchSets{cur});
        
        l_objectReflectances = [l_objectReflectances; l_curObjectReflectances(2 : end, :)];
        
    end
    
else 
    
    l_objectReflectances = IDTTaskData.IDTCreationConstraints_In_PatchSet;
    
end

HDM_OFT_Utils.OFT_DispTitle('start IDT profile creation');

for curIlluminandIndex=1:size(OFT_Illuminations,2)
    
    TargetIlluminantStr=OFT_Illuminations{curIlluminandIndex};
    HDM_OFT_Utils.OFT_DispSubTitle(strcat('create idt profile for illuminant: ', TargetIlluminantStr));
    
    [OFT_IDTFile, OFT_IDT_B, OFT_IDT_b]=...
        HDM_OFT_IDT_MinimumError...
            (l_objectReflectances,...
            IDTTaskData.SpectralResponse_Out_SpectralResponseFile,...
            IDTTaskData.PreLinearisation_Out_LinCurve,...
            HDM_OFT_ColorNeutralCompensations.ChromaticAdaptationBradfordType(),...
            IDTTaskData.IDTCreationConstraints_In_ErrorMinimizationDomain,...
            TargetIlluminantStr, IDTTaskData.IDTCreationConstraints_In_SceneIllumination,...
            ...% HDM_OFT_IDT_ReferenceCamera.CIEType(),HDM_OFT_CIEStandard.StandardObserver1931_2Degrees());
            HDM_OFT_IDT_ReferenceCamera.RICDType(),HDM_OFT_CIEStandard.StandardObserver1931_2Degrees());

    OFT_IDTFiles{curIlluminandIndex}=OFT_IDTFile;    

    disp(strcat('idt file: ', OFT_IDTFile));
    disp('estimated matrix');
    disp(OFT_IDT_B);
    
    global gOFT_SpectralDataBasedTransformedImage2View;
    
    if(isempty(gOFT_SpectralDataBasedTransformedImage2View))
    
        gOFT_SpectralDataBasedTransformedImage2View = HDM_OFT_EvaluateIDTProfiledChartImage...
        (OFT_IDT_B, OFT_IDT_b, ...
        IDTTaskData.Evaluation_In_TestImage, IDTTaskData.PreLinearisation_Out_LinCurve, ...
        TargetIlluminantStr, IDTTaskData.IDTCreationConstraints_In_SceneIllumination, ...
        HDM_OFT_PatchSet.GretagMacbethColorChecker(),... 
        IDTTaskData.SpectralResponse_Out_SpectralResponseFile, ...
        HDM_OFT_IDT_ReferenceCamera.RICDType(), HDM_OFT_ColorNeutralCompensations.ChromaticAdaptationBradfordType(), HDM_OFT_CIEStandard.StandardObserver1931_2Degrees());

    else

        gOFT_SpectralDataBasedTransformedImage2View = HDM_OFT_EvaluateIDTProfiledChartImage...
        (OFT_IDT_B, OFT_IDT_b, ...
        IDTTaskData.Evaluation_In_TestImage, IDTTaskData.PreLinearisation_Out_LinCurve, ...
        TargetIlluminantStr, IDTTaskData.IDTCreationConstraints_In_SceneIllumination, ...
        HDM_OFT_PatchSet.GretagMacbethColorChecker(),... 
        IDTTaskData.SpectralResponse_Out_SpectralResponseFile, ...
        HDM_OFT_IDT_ReferenceCamera.RICDType(), HDM_OFT_ColorNeutralCompensations.ChromaticAdaptationBradfordType(), HDM_OFT_CIEStandard.StandardObserver1931_2Degrees(), ...
        gOFT_SpectralDataBasedTransformedImage2View);

    end
    
end

HDM_OFT_Utils.OFT_DispTitle('finish IDT profile creation');

end