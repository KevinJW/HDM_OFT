function [o_B, o_resnormBEstimation] = HDM_OFT_Estimate3x3TransformationMatrixConstrained ...
    (i_ErrorMinimizationDomain, ...
    i_ReferenceTristimuli, i_CurrentTristimuli, i_M, i_w, ...
    i_PreserveColourPositions, i_weightingVec)

    %%todo n equalities

    o_B = [];
    o_resnormBEstimation = [];

    l_B0 = ...
        [1 0 0;
        0 1 0;
        0 0 1];
    
    % transpose as is usual in linear algebra:
    l_CurrentTristimuli = i_CurrentTristimuli';
    l_ReferenceTristimuli = i_ReferenceTristimuli';
    
    l_refWhite2Preserve = l_ReferenceTristimuli(i_PreserveColourPositions, :);
    l_curWhite2Preserve = l_CurrentTristimuli(i_PreserveColourPositions, :);
    
    l_curWhite = [];
    
    for cur = 1 : size(l_curWhite2Preserve, 1)
    
        l_a = l_curWhite2Preserve(cur, 1);
        l_b = l_curWhite2Preserve(cur, 2);
        l_c = l_curWhite2Preserve(cur, 3);


        l_cur = [ l_a, l_b, l_c, 0, 0, 0, 0, 0, 0; ...
                    0, 0, 0, l_a, l_b, l_c, 0, 0, 0; ...
                    0, 0, 0, 0, 0, 0, l_a, l_b, l_c   ];

        l_curWhite = [l_curWhite; l_cur];
            
    end       
    
    l_W = diag(sqrt(i_weightingVec));
    l_CurrentWeigthedTristimuli = l_W * l_CurrentTristimuli;
    l_ReferenceWeigthedTristimuli = l_W * l_ReferenceTristimuli;

    switch i_ErrorMinimizationDomain
        case 'Lab'  
            
            l_ReferenceTristimuli_LabDomain = HDM_OFT_ColorConversions.OFT_CIELab(i_ReferenceTristimuli, i_w);           
            
            if isempty(i_weightingVec)
                
                [l_x, l_res] = fmincon(@(l_B0) LSFunConLab(l_B0, l_ReferenceTristimuli_LabDomain, i_CurrentTristimuli, i_w, i_M), ...
                                l_B0, ...
                                [], [], ...
                                l_curWhite, ...
                                l_refWhite2Preserve', ...
                                [],[]);
                            
            else

                l_weightedRef = l_W * l_ReferenceTristimuli_LabDomain';
                
                [l_x, l_res] = fmincon(@(l_B0) LSFunConLab_Weighted(l_B0, l_weightedRef, i_CurrentTristimuli, i_w, i_M, l_W), ...
                                l_B0, ...
                                [], [], ...
                                l_curWhite, ...
                                l_refWhite2Preserve', ...
                                [],[]);
                            
                % for check
                
            end
                        
            o_B = l_x;
            o_resnormBEstimation = l_res;
                        
            %for check
            l_refWhite2PreserveCalculated = l_curWhite2Preserve * l_x;                       
            

        case 'Luv'
            
            l_ReferenceTristimuli_LuvDomain = HDM_OFT_ColorConversions.OFT_CIELuv(i_ReferenceTristimuli, i_w);           

            if isempty(i_weightingVec)
                
                [l_x, l_res] = fmincon(@(l_B0) LSFunConLuv(l_B0, l_ReferenceTristimuli_LuvDomain, i_CurrentTristimuli, i_w, i_M), ...
                                l_B0, ...
                                [], [], ...
                                l_curWhite, ...
                                l_refWhite2Preserve', ...
                                [],[]);
                        
            else
                
                l_weightedRef = l_W * l_ReferenceTristimuli_LuvDomain';
                
                [l_x, l_res] = fmincon(@(l_B0) LSFunConLuv_Weighted(l_B0, l_weightedRef, i_CurrentTristimuli, i_w, i_M, l_W), ...
                                l_B0, ...
                                [], [], ...
                                l_curWhite, ...
                                l_refWhite2Preserve', ...
                                [],[]);
                
            end
                                                          
            o_B = l_x;
            o_resnormBEstimation = l_res;
            
            %for check
            l_refWhite2PreserveCalculated = l_curWhite2Preserve * l_x;           
            
        case 'XYZ'    
            
            if isempty(i_weightingVec)
                
                [l_WPPLSresSep, l_WPPLSresSep_resNorm] = OFT_lsqlinNDim_WithConstraints ...
                    (l_B0, l_ReferenceTristimuli, l_CurrentTristimuli, l_refWhite2Preserve, l_curWhite2Preserve);
            
            else
                
                [l_WPPLSresSep, l_WPPLSresSep_resNorm] = OFT_lsqlinNDim_WithConstraints ...
                    (l_B0, l_ReferenceWeigthedTristimuli, l_CurrentWeigthedTristimuli, l_refWhite2Preserve, l_curWhite2Preserve);
            
            end
            
            o_B = l_WPPLSresSep;
            o_resnormBEstimation = l_WPPLSresSep_resNorm;
            
            % for check
            l_uncon = (l_CurrentTristimuli' * l_CurrentTristimuli)^(-1) * l_CurrentTristimuli' * l_ReferenceTristimuli;           
            l_wComputed = l_curWhite2Preserve * l_WPPLSresSep;
            
        otherwise
    end
    
    o_B(o_B < 0.0001) = 0.0;
    
    o_B = o_B';

end

%%CAM Models

function o_out = LSFunConLab(i_B, i_refVals, i_srcVals, i_w, i_M)

    l_D = ...
        i_refVals - ...
        HDM_OFT_ColorConversions.OFT_CIELab(i_M * i_B * i_srcVals, i_w);

    l_norm = sqrt(sum(abs(l_D).^2, 1)); % or delta E 2000

    o_out = l_norm';
    
    o_out = sqrt(sum(l_norm));

end

function o_out = LSFunConLab_Weighted(i_B, i_refVals, i_srcVals, i_w, i_M, i_W)

    l_D = ...
        i_refVals - ...
        (i_W * HDM_OFT_ColorConversions.OFT_CIELab(i_M * i_B * i_srcVals, i_w)');

    l_norm = sqrt(sum(abs(l_D).^2, 1)); % or delta E 2000

    o_out = l_norm';
    
    o_out = sqrt(sum(l_norm));

end

function o_out = LSFunConLuv(i_B, i_refVals, i_srcVals, i_w, i_M)

    l_D = ...
        i_refVals - ...
        HDM_OFT_ColorConversions.OFT_CIELuv(i_M * i_B * i_srcVals, i_w);

    l_norm = sqrt(sum(abs(l_D).^2, 1));

    o_out = l_norm';
    
    o_out = sqrt(sum(l_norm));

end

function o_out = LSFunConLuv_Weighted(i_B, i_refVals, i_srcVals, i_w, i_M, i_W)

    l_D = ...
        i_refVals - ...
        (i_W * HDM_OFT_ColorConversions.OFT_CIELuv(i_M * i_B * i_srcVals, i_w)');

    l_norm = sqrt(sum(abs(l_D).^2, 1));

    o_out = l_norm';
    
    o_out = sqrt(sum(l_norm));

end

%% Primary Models

%% utilities

function [o_WPPLSresSep, o_WPPLSresSep_resNorm] = OFT_lsqlinNDim_WithConstraints ...
    (i_B0, i_ReferenceTristimuli, i_CurrentTristimuli, i_refWhite2Preserve, i_curWhite2Preserve)

    %% solve constraint linear LS for each vector independently, see above
    [l_WPPLSresR, l_WPPLSresR_resNorm] = lsqlin(i_CurrentTristimuli, i_ReferenceTristimuli(: , 1), ...
        [], [], i_curWhite2Preserve, i_refWhite2Preserve(:, 1), [],[]);

    [l_WPPLSresG, l_WPPLSresG_resNorm] = lsqlin(i_CurrentTristimuli, i_ReferenceTristimuli(: , 2), ...
        [], [], i_curWhite2Preserve, i_refWhite2Preserve(:, 2), [],[]);

    [l_WPPLSresB, l_WPPLSresB_resNorm] = lsqlin(i_CurrentTristimuli, i_ReferenceTristimuli(: , 3), ...
        [], [], i_curWhite2Preserve, i_refWhite2Preserve(:, 3), [],[]);

    %% combine results
    o_WPPLSresSep = [l_WPPLSresR, l_WPPLSresG, l_WPPLSresB];
    o_WPPLSresSep_resNorm = [l_WPPLSresR_resNorm, l_WPPLSresG_resNorm, l_WPPLSresB_resNorm];

end

