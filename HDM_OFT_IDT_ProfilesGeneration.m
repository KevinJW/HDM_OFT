function OFT_IDTFiles=HDM_OFT_IDT_ProfilesGeneration...
    (i_IDTTaskData)

if isempty(i_IDTTaskData.SpectralResponse_In_LightCalibrationSpectrum)        %chart based
    OFT_Illuminations={i_IDTTaskData.IDTCreationConstraints_In_WhitePoint};
else                                                                        %spectral based therefore we use the line and light spectra too
    OFT_Illuminations={i_IDTTaskData.IDTCreationConstraints_In_WhitePoint};%, ...
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

if not(isempty(i_IDTTaskData.IDTCreationConstraints_In_AdditionalPatchSets))
    
    l_objectReflectances = HDM_OFT_PatchSet.GetPatchSpectra(i_IDTTaskData.IDTCreationConstraints_In_PatchSet);
    
    if not(strcmp(i_IDTTaskData.IDTCreationConstraints_In_PatchSet_Illuminant, 'E'))
        
        l_illuminantSpectrum4PatchSet = HDM_OFT_GetIlluminantSpectrum(i_IDTTaskData.IDTCreationConstraints_In_PatchSet_Illuminant);

        for cur = 2 : size(l_objectReflectances, 1)

            l_objectReflectances(cur) = l_objectReflectances(cur) ./ l_illuminantSpectrum4PatchSet(2);

        end

    end
    
    for cur = 1 : size(i_IDTTaskData.IDTCreationConstraints_In_AdditionalPatchSets, 1)
        
        l_curObjectReflectances = HDM_OFT_PatchSet.GetPatchSpectra(i_IDTTaskData.IDTCreationConstraints_In_AdditionalPatchSets{cur});
        
        l_illuminatStr = i_IDTTaskData.IDTCreationConstraints_In_AdditionalPatchSets_Illuminant{cur};
        if not(strcmp(l_illuminatStr, 'E'))
            
            l_illuminantSpectrum4PatchSet = HDM_OFT_GetIlluminantSpectrum(l_illuminatStr);
            
            for icur = 2 : size(l_curObjectReflectances, 1)
                
                l_curObjectReflectances(icur) = l_curObjectReflectances(icur) ./ l_illuminantSpectrum4PatchSet(2);
                               
            end
                                  
        end
        
        l_curObjectReflectances(isnan(l_curObjectReflectances)) = 0; 
        
        l_objectReflectances = [l_objectReflectances; l_curObjectReflectances(2 : end, :)];
        
    end
    
else 
    
    l_objectReflectances = i_IDTTaskData.IDTCreationConstraints_In_PatchSet;
    
    if not(strcmp(i_IDTTaskData.IDTCreationConstraints_In_PatchSet_Illuminant, 'E'))
        
        l_objectReflectances = HDM_OFT_PatchSet.GetPatchSpectra(i_IDTTaskData.IDTCreationConstraints_In_PatchSet);
            
        l_illuminantSpectrum4PatchSet = HDM_OFT_GetIlluminantSpectrum(i_IDTTaskData.IDTCreationConstraints_In_PatchSet_Illuminant);

        for cur = 2 : size(l_objectReflectances, 1)

            l_objectReflectances(cur) = l_objectReflectances(cur) ./ l_illuminantSpectrum4PatchSet(2);

        end

    end
    
end

HDM_OFT_Utils.OFT_DispTitle('start IDT profile creation');

for curIlluminandIndex=1:size(OFT_Illuminations, 2)
    
    l_targetIlluminantStr=OFT_Illuminations{curIlluminandIndex};
    HDM_OFT_Utils.OFT_DispSubTitle(strcat('create idt profile for illuminant: ', l_targetIlluminantStr));
    
    l_neutralsCompensation = HDM_OFT_ColorNeutralCompensations.ChromaticAdaptationBradfordType();
    l_errorMinimizationDomain = i_IDTTaskData.IDTCreationConstraints_In_ErrorMinimizationDomain;
    l_referenceDomain = HDM_OFT_IDT_ReferenceCamera.RICDType();
    
    [l_IDT_B, l_IDT_b, l_referenceWhite, l_resnormBEstimation, l_meanDeltaE2000, l_stdDevDeltaE2000]=...
        HDM_OFT_IDT_MinimumError...
            (l_objectReflectances,...
            i_IDTTaskData.SpectralResponse_Out_SpectralResponseFile,...
            i_IDTTaskData.PreLinearisation_Out_LinCurve,...
            l_neutralsCompensation,...
            l_errorMinimizationDomain,...
            l_targetIlluminantStr, i_IDTTaskData.IDTCreationConstraints_In_SceneIllumination,...
            ...% HDM_OFT_IDT_ReferenceCamera.CIEType(),HDM_OFT_CIEStandard.StandardObserver1931_2Degrees());
            l_referenceDomain, HDM_OFT_CIEStandard.StandardObserver1931_2Degrees());
        
        
    %% write idt profile file and append results to statistics
    HDM_OFT_Utils.OFT_DispTitle('write idt and icc profile file and append results to statistics');
    [l_IDT_File, l_ICC_File] = HDM_OFT_WriteProfileAndStatEntry...
    (i_IDTTaskData.SpectralResponse_Out_SpectralResponseFile, l_resnormBEstimation, l_IDT_B, l_IDT_b,...
    l_neutralsCompensation, l_targetIlluminantStr, l_referenceWhite, l_referenceDomain, l_errorMinimizationDomain, ...
    l_meanDeltaE2000, l_stdDevDeltaE2000, i_IDTTaskData);


    OFT_IDTFiles{(2 * (curIlluminandIndex - 1)) + 1} = l_IDT_File;    
    OFT_IDTFiles{(2 * (curIlluminandIndex - 1))+ 2} = l_ICC_File;

    disp(strcat('idt file: ', l_IDT_File));
    disp('estimated matrix');
    disp(l_IDT_B);
    
    global gOFT_SpectralDataBasedTransformedImage2View;
    
    if(isempty(gOFT_SpectralDataBasedTransformedImage2View))
    
        gOFT_SpectralDataBasedTransformedImage2View = HDM_OFT_EvaluateIDTProfiledChartImage...
        (l_IDT_B, l_IDT_b, ...
        i_IDTTaskData.Evaluation_In_TestImage, i_IDTTaskData.PreLinearisation_Out_LinCurve, ...
        l_targetIlluminantStr, i_IDTTaskData.IDTCreationConstraints_In_SceneIllumination, ...
        HDM_OFT_PatchSet.GretagMacbethColorChecker(),... 
        i_IDTTaskData.SpectralResponse_Out_SpectralResponseFile, ...
        HDM_OFT_IDT_ReferenceCamera.RICDType(), HDM_OFT_ColorNeutralCompensations.ChromaticAdaptationBradfordType(), HDM_OFT_CIEStandard.StandardObserver1931_2Degrees());

    else

        gOFT_SpectralDataBasedTransformedImage2View = HDM_OFT_EvaluateIDTProfiledChartImage...
        (l_IDT_B, l_IDT_b, ...
        i_IDTTaskData.Evaluation_In_TestImage, i_IDTTaskData.PreLinearisation_Out_LinCurve, ...
        l_targetIlluminantStr, i_IDTTaskData.IDTCreationConstraints_In_SceneIllumination, ...
        HDM_OFT_PatchSet.GretagMacbethColorChecker(),... 
        i_IDTTaskData.SpectralResponse_Out_SpectralResponseFile, ...
        HDM_OFT_IDT_ReferenceCamera.RICDType(), HDM_OFT_ColorNeutralCompensations.ChromaticAdaptationBradfordType(), HDM_OFT_CIEStandard.StandardObserver1931_2Degrees(), ...
        gOFT_SpectralDataBasedTransformedImage2View);

    end
    
end

HDM_OFT_Utils.OFT_DispTitle('finish IDT profile creation');

end