function [o_B, o_resnormBEstimation] = HDM_OFT_Estimate3x3TransformationMatrixConstrained ...
    (i_ErrorMinimizationDomain, ...
    i_ReferenceTristimuli, i_CurrentTristimuli, i_M, i_w, ...
    i_weightingVec, i_PreserveColourPositions)

    %%todo n equalities

    o_B = [];
    o_resnormBEstimation = [];

    l_B0 = ...
        [1 0;
        0 1;
        0 0];
    
    l_refWhite2Preserve = i_ReferenceTristimuli(i_PreserveColourPositions, :);
    
    l_curWhite2Preserve = i_CurrentTristimuli(i_PreserveColourPositions, :);
    
    l_curWhite = [];
    
    for cur = 1 : size(l_curWhite2Preserve, 1)
    
        l_a = l_curWhite2Preserve(cur, 1);
        l_b = l_curWhite2Preserve(cur, 2);
        l_c = l_curWhite2Preserve(cur, 3);

        l_cur = [ l_a, l_b, l_c, 0, 0, 0, 0, 0, 0; ...
                    0, 0, 0, l_a, l_b, l_c, 0, 0, 0; ...
                    0, 0, 0, 0, 0, 0, l_a, l_b, l_c   ];
                
        l_cur = [ l_a, l_b, 0, 0, 0, 0; ...
                    0, 0, l_a, l_b, 0, 0; ...
                    0, 0, 0, 0, l_a, l_b];

        l_curWhite = [l_curWhite; l_cur];
            
    end       
    
    if not(exist('i_weightingVec', 'var'))

        i_weightingVec = [];
        
    end
    
    if isempty(i_weightingVec)

        l_ones = zeros(size(i_CurrentTristimuli,1), 1);
        l_ones(:, 1) = 1;
        l_W = diag(l_ones);
        
    else
        
        l_W = diag(sqrt(i_weightingVec));
        
    end
    
    l_CurrentWeigthedTristimuli = l_W * i_CurrentTristimuli;
    l_ReferenceWeigthedTristimuli = l_W * i_ReferenceTristimuli;

    switch i_ErrorMinimizationDomain
        case 'Lab'  
            
            l_ReferenceTristimuli_LabDomain = HDM_OFT_ColorConversions.OFT_CIELab(i_ReferenceTristimuli', i_w)';           
            
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
            
            l_B = zeros(3, 3);
            
            l_B(1, 1) = o_B(1, 1);
            l_B(2, 1) = o_B(2, 1);
            l_B(3, 1) = o_B(3, 1);

            l_B(1, 2) = o_B(1, 2);
            l_B(2, 2) = o_B(2, 2);
            l_B(3, 2) = o_B(3, 2);

            l_B(1, 3) = 1 - o_B(1, 1) - o_B(1, 2);
            l_B(2, 3) = 1 - o_B(2, 1) - o_B(2, 2);
            l_B(3, 3) = 1 - o_B(3, 1) - o_B(3, 2);
            
            o_B = l_B;                       
            
            l_refWhite2Preserve;
            l_MB = i_M * o_B;
            l_wT = l_curWhite2Preserve * o_B;
            l_wT = l_curWhite2Preserve * l_MB;
                        
            %for check
            l_refWhite2PreserveCalculated = l_curWhite2Preserve * l_x;                       
            

        case 'Luv'
            
            l_ReferenceTristimuli_LuvDomain = HDM_OFT_ColorConversions.OFT_CIELuv(i_ReferenceTristimuli', i_w)';           

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
                    (l_B0, i_ReferenceTristimuli, i_CurrentTristimuli, l_refWhite2Preserve, l_curWhite2Preserve);
            
            else
                
                [l_WPPLSresSep, l_WPPLSresSep_resNorm] = OFT_lsqlinNDim_WithConstraints ...
                    (l_B0, l_ReferenceWeigthedTristimuli, l_CurrentWeigthedTristimuli, l_refWhite2Preserve, l_curWhite2Preserve);
            
            end
            
            o_B = l_WPPLSresSep;
            o_resnormBEstimation = l_WPPLSresSep_resNorm;
            
            l_refWhite2Preserve;
            l_MB = i_M * o_B;
            l_wT = l_curWhite2Preserve * o_B;
            %l_wT = l_MB * l_curWhite2Preserve';
            
            % for check
            l_uncon = (i_CurrentTristimuli' * i_CurrentTristimuli)^(-1) * i_CurrentTristimuli' * i_ReferenceTristimuli;           
            l_wComputed = l_curWhite2Preserve * l_WPPLSresSep;
            
        otherwise
    end
    
    o_B(abs(o_B) < 0.0001) = 0.0;

end

%%CAM Models

function o_out = LSFunConLab(i_B, i_refVals, i_srcVals, i_w, i_M)

    l_rangeLB = 1;
    l_rangeUB = 24;
    
    l_sampleSelection = [l_rangeLB : l_rangeUB];
    
    l_B = i_B;
    
    l_B(1, 1) = i_B(1, 1);
    l_B(2, 1) = i_B(2, 1);
    l_B(3, 1) = i_B(3, 1);
    
    l_B(1, 2) = i_B(1, 2);
    l_B(2, 2) = i_B(2, 2);
    l_B(3, 2) = i_B(3, 2);
    
    l_B(1, 3) = 1 - i_B(1, 1) - i_B(1, 2);
    l_B(2, 3) = 1 - i_B(2, 1) - i_B(2, 2);
    l_B(3, 3) = 1 - i_B(3, 1) - i_B(3, 2);
    
    l_D = ...
        i_refVals(l_sampleSelection, :) - ...
        HDM_OFT_ColorConversions.OFT_CIELab(i_M * l_B * i_srcVals(l_sampleSelection, :)', i_w)';

    l_norm = sqrt(sum(abs(l_D).^2, 1));

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

    l_lbr = [];%[0.3 0.05 0];
    l_ubr = [];%[4 4 4];
    
    l_rangeLB = 1;
    l_rangeUB = 24;
    
    l_sampleSelection = [l_rangeLB : l_rangeUB];
    
    %l_sampleSelection = [1 : 5];

    %% solve constraint linear LS for each vector independently, see above
    [l_WPPLSresR, l_WPPLSresR_resNorm] = lsqlin(i_CurrentTristimuli(l_sampleSelection, :), i_ReferenceTristimuli(l_sampleSelection, 1), ...
        [], [], i_curWhite2Preserve, i_refWhite2Preserve(:, 1), l_lbr, l_ubr);
    
    l_lbg = [];%[0.05 0.1 0]
    l_ubg = [];%[4 4 4];

    [l_WPPLSresG, l_WPPLSresG_resNorm] = lsqlin(i_CurrentTristimuli(l_sampleSelection, :), i_ReferenceTristimuli(l_sampleSelection, 2), ...
        [], [], i_curWhite2Preserve, i_refWhite2Preserve(:, 2), l_lbg, l_ubg);
    
    l_lbb = [];%[0.05 0.05 0]
    l_ubb = [];%[0.3 0.3 4];

    [l_WPPLSresB, l_WPPLSresB_resNorm] = lsqlin(i_CurrentTristimuli(l_sampleSelection, :), i_ReferenceTristimuli(l_sampleSelection, 3), ...
        [], [], i_curWhite2Preserve, i_refWhite2Preserve(:, 3), l_lbb, l_ubb);

    %% combine results
    o_WPPLSresSep = [l_WPPLSresR, l_WPPLSresG, l_WPPLSresB];
    o_WPPLSresSep_resNorm = [l_WPPLSresR_resNorm, l_WPPLSresG_resNorm, l_WPPLSresB_resNorm];

end

