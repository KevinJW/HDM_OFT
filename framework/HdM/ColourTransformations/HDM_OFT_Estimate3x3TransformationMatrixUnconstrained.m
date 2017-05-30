function [o_B, o_resnormBEstimation] = HDM_OFT_Estimate3x3TransformationMatrixUnconstrained ...
    (i_ErrorMinimizationDomain, ...
    i_ReferenceTristimuli, i_CurrentTristimuli, i_M, i_w, i_weightingVec)

    l_B0 = ...
        [1 0 0;
        0 1 0;
        0 0 1];
    
    % transpose as is usual in linear algebra:
    l_CurrentTristimuli = i_CurrentTristimuli';
    l_ReferenceTristimuli = i_ReferenceTristimuli';
    
    l_M = i_M';
    
    if not(isempty(i_weightingVec))
        
        l_W = diag(sqrt(i_weightingVec));
        l_CurrentWeigthedTristimuli = l_W * l_CurrentTristimuli;
        l_ReferenceWeigthedTristimuli = l_W * l_ReferenceTristimuli;
        
    end

    switch i_ErrorMinimizationDomain
        case 'Lab'           
            
            l_ReferenceTristimuli_LabDomain = HDM_OFT_ColorConversions.OFT_CIELab(i_ReferenceTristimuli, i_w);
            
            if isempty(i_weightingVec)
                
                [o_B,o_resnormBEstimation] = ...
                    lsqnonlin( ...
                    @(l_B0) OFT_IDT_MeritFunctionCoreLab(l_B0, l_ReferenceTristimuli_LabDomain, i_CurrentTristimuli, i_M, i_w), ...
                    l_B0);
            else
                
                l_weightedRef = l_W * l_ReferenceTristimuli_LabDomain';

                [o_B,o_resnormBEstimation] = ...
                    lsqnonlin( ...
                    @(l_B0) OFT_IDT_MeritFunctionCoreLab_Weighted(l_B0, l_weightedRef', i_CurrentTristimuli, i_M, i_w,  l_W), ...
                    l_B0);
            
            end
            
        case 'Luv'
            
            l_ReferenceTristimuli_LuvDomain = HDM_OFT_ColorConversions.OFT_CIELuv(i_ReferenceTristimuli, i_w);
            
            if isempty(i_weightingVec)
                
                [o_B,o_resnormBEstimation] = ...
                    lsqnonlin( ...
                    @(l_B0) OFT_IDT_MeritFunctionCoreLuv(l_B0, l_ReferenceTristimuli_LuvDomain, i_CurrentTristimuli, i_M, i_w), ...
                    l_B0);
            
            else
                
                l_weightedRef = l_W * l_ReferenceTristimuli_LuvDomain';

                [o_B,o_resnormBEstimation] = ...
                    lsqnonlin( ...
                    @(l_B0) OFT_IDT_MeritFunctionCoreLuv_Weighted(l_B0, l_weightedRef', i_CurrentTristimuli, i_M, i_w,  l_W), ...
                    l_B0);
            
            end
            
        case 'XYZ'                               
            
            % base matlab lin problem: pseudo inverse, ls if over determined
            % M = N x or M^T = x^T N^T
            
            if isempty(i_weightingVec)
                
                l_B1 = (l_CurrentTristimuli' * l_CurrentTristimuli)^(-1) * l_CurrentTristimuli' * l_ReferenceTristimuli;
            
            else
                
                l_B1 = (l_CurrentWeigthedTristimuli' * l_CurrentWeigthedTristimuli)^(-1) * ...
                    l_CurrentWeigthedTristimuli' * l_ReferenceWeigthedTristimuli;
            
            end
                                   
            o_B = l_B1;
            o_resnormBEstimation = [];
            
            %% test other impl

            % base matlab lin problem: recommended, ls if over determined
            l_B2 = mldivide(l_CurrentTristimuli, l_ReferenceTristimuli);

            % base matlab lin problem: alternative, ls if over determined
            l_B3 = linsolve(l_CurrentTristimuli, l_ReferenceTristimuli);
            
            %lsqlin with constraints must be used separately
            l_B3a1 = lsqlin(l_CurrentTristimuli, l_ReferenceTristimuli(: , 1), ...
                [], [], [], [], [],[]);
            l_B3a2 = lsqlin(l_CurrentTristimuli, l_ReferenceTristimuli(: , 2), ...
                [], [], [], [], [],[]);
            l_B3a3 = lsqlin(l_CurrentTristimuli, l_ReferenceTristimuli(: , 3), ...
                [], [], [], [], [],[]);
            
            l_B3a = [l_B3a1, l_B3a2, l_B3a3];
            
            % optim toolbox: can be used also, although it is a linear problem
            if isempty(i_weightingVec)
                
                [l_B4, o_resnormBEstimation] = ...
                    lsqnonlin(@(l_B0) OFT_IDT_MeritFunctionCoreXYZ(l_B0, l_ReferenceTristimuli, l_CurrentTristimuli), l_B0);
            
            else
                
                [l_B4W, o_resnormBEstimationW] = ...
                    lsqnonlin(@(l_B0) OFT_IDT_MeritFunctionCoreXYZ_Weighted ...
                    (l_B0, l_ReferenceWeigthedTristimuli, l_CurrentTristimuli, l_W), l_B0);
            
            end
            
        otherwise
    end
    
    o_B(o_B < 0.0001) = 0.0;
    
    o_B = o_B';

end

%%CAM Models

function o_F = OFT_IDT_MeritFunctionCoreLuv(i_B0, i_ReferenceTristimuli, i_CurrentTristimuli, i_M, i_w)

k = 1 : size(i_ReferenceTristimuli,2);

o_F = i_ReferenceTristimuli(:,k) -...
        HDM_OFT_ColorConversions.OFT_CIELuv(i_M * i_B0 * i_CurrentTristimuli(:,k), i_w);

l_norm = sqrt(sum(abs(o_F).^2, 1));

o_F = l_norm';

end

function o_F = OFT_IDT_MeritFunctionCoreLuv_Weighted(i_B0, i_ReferenceTristimuli, i_CurrentTristimuli, i_M, i_w, i_W)

k = 1 : size(i_ReferenceTristimuli,2);

o_F = i_ReferenceTristimuli(:,k) -...
        (i_W * HDM_OFT_ColorConversions.OFT_CIELuv(i_M * i_B0 * i_CurrentTristimuli(:,k), i_w)')';

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

function o_F = OFT_IDT_MeritFunctionCoreLab_Weighted(i_B0, i_ReferenceTristimuli, i_CurrentTristimuli, i_M, i_w, i_W)

k = 1 : size(i_ReferenceTristimuli,2);

o_F = i_ReferenceTristimuli(:,k) - ...
        (i_W * HDM_OFT_ColorConversions.OFT_CIELab(i_M * i_B0 * i_CurrentTristimuli(:,k), i_w)')';

l_norm = sqrt(sum(abs(o_F).^2, 1)); % or delta E 2000

o_F = l_norm';

end

%% Primary Models

function o_F = OFT_IDT_MeritFunctionCoreXYZ(i_B0, i_ReferenceTristimuli, i_CurrentTristimuli)

o_F = i_ReferenceTristimuli - i_CurrentTristimuli * i_B0;

end

function o_F = OFT_IDT_MeritFunctionCoreXYZ_Weighted(i_B0, i_ReferenceTristimuliWeighted, i_CurrentTristimuli, i_W)

o_F = i_ReferenceTristimuliWeighted - i_W * i_CurrentTristimuli * i_B0;

end

