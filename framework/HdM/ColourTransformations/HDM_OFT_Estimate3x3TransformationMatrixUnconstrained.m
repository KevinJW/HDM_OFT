function [o_B, o_resnormBEstimation] = HDM_OFT_Estimate3x3TransformationMatrixUnconstrained ...
    (i_ErrorMinimizationDomain, ...
    i_ReferenceTristimuli, i_CurrentTristimuli, i_M, i_w)

    l_B0 = ...
        [1 0 0;
        0 1 0;
        0 0 1];
    
    l_k = 1 : size(i_ReferenceTristimuli,2);

    switch i_ErrorMinimizationDomain
        case 'Lab'           
            
            l_ReferenceTristimuli_LabDomain = HDM_OFT_ColorConversions.OFT_CIELab(i_ReferenceTristimuli(:, l_k), i_w);
            
            [o_B,o_resnormBEstimation] = ...
                lsqnonlin(@(l_B0) OFT_IDT_MeritFunctionCoreLab(l_B0, l_ReferenceTristimuli_LabDomain, i_CurrentTristimuli, i_M, i_w), l_B0);
            
        case 'Luv'
            
            l_ReferenceTristimuli_LuvDomain = HDM_OFT_ColorConversions.OFT_CIELuv(i_ReferenceTristimuli(:, l_k), i_w);
            
            [o_B,o_resnormBEstimation] = ...
                lsqnonlin(@(l_B0) OFT_IDT_MeritFunctionCoreLuv(l_B0, l_ReferenceTristimuli_LuvDomain, i_CurrentTristimuli, i_M, i_w), l_B0);
            
        case 'XYZ'
            
            [o_B,o_resnormBEstimation] = ...
                lsqnonlin(@(l_B0) OFT_IDT_MeritFunctionCoreXYZ(l_B0, i_ReferenceTristimuli, i_CurrentTristimuli), l_B0);
            
        otherwise
    end

end

function o_F = OFT_IDT_MeritFunctionCoreLuv(i_B0, i_ReferenceTristimuli, i_CurrentTristimuli, i_M, i_w)

k = 1 : size(i_ReferenceTristimuli,2);

o_F = i_ReferenceTristimuli(:,k) -...
        HDM_OFT_ColorConversions.OFT_CIELuv(i_M * i_B0 * i_CurrentTristimuli(:,k), i_w);

l_norm = sqrt(sum(abs(o_F).^2, 1));

o_F = l_norm';

end

function o_F = OFT_IDT_MeritFunctionCoreLab(i_B0, i_ReferenceTristimuli, i_CurrentTristimuli, i_M, i_w)

k = 1 : size(i_ReferenceTristimuli,2);

o_F = i_ReferenceTristimuli(:,k) - ...
        HDM_OFT_ColorConversions.OFT_CIELab(i_M * i_B0 * i_CurrentTristimuli(:,k), i_w);

l_norm = sqrt(sum(abs(o_F).^2, 1)); % or delta E 2000

o_F = l_norm';

end

function o_F = OFT_IDT_MeritFunctionCoreXYZ(i_B0, i_ReferenceTristimuli, i_CurrentTristimuli)

k = 1 : size(i_ReferenceTristimuli,2);

o_F = i_ReferenceTristimuli(:,k) - i_B0 * i_CurrentTristimuli(:,k);

l_norm = sqrt(sum(abs(o_F).^2, 1));

o_F = l_norm';

end

