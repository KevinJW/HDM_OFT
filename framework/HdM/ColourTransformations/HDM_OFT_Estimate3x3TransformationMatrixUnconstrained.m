function [o_B, o_resnormBEstimation] = HDM_OFT_Estimate3x3TransformationMatrixUnconstrained ...
    (i_ErrorMinimizationDomain, ...
    i_ReferenceTristimuli, i_CurrentTristimuli, i_M, i_w, i_weightingVec)

    l_B0 = ...
        [1 0;
        0 1;
        0 0];
    
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
    
    options = optimoptions('lsqnonlin');
    options.MaxFunEvals = 4 * 600;
    

    switch i_ErrorMinimizationDomain
        case 'Lab'           
            
            l_ReferenceTristimuli_LabDomain = HDM_OFT_ColorConversions.OFT_CIELab(i_ReferenceTristimuli', i_w)';
            
            if isempty(i_weightingVec)
                
                [o_B,o_resnormBEstimation] = ...
                    lsqnonlin( ...
                    @(l_B0) OFT_IDT_MeritFunctionCoreLab(l_B0, l_ReferenceTristimuli_LabDomain, i_CurrentTristimuli, i_M, i_w), ...
                    l_B0, [], [], options);
            else
                
                l_weightedRef = l_W * l_ReferenceTristimuli_LabDomain';

                [o_B,o_resnormBEstimation] = ...
                    lsqnonlin( ...
                    @(l_B0) OFT_IDT_MeritFunctionCoreLab_Weighted(l_B0, l_weightedRef', i_CurrentTristimuli, i_M, i_w,  l_W), ...
                    l_B0);
            
            end
            
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
            
        case 'Luv'
            
            l_ReferenceTristimuli_LuvDomain = HDM_OFT_ColorConversions.OFT_CIELuv(i_ReferenceTristimuli', i_w)';
            
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
            
        case 'XYZ'                               
            
            % base matlab lin problem: pseudo inverse, ls if over determined
            % M = N x or M^T = x^T N^T
            
            if isempty(i_weightingVec)
                
                l_B1 = (i_CurrentTristimuli' * i_CurrentTristimuli)^(-1) * i_CurrentTristimuli' * i_ReferenceTristimuli;
            
            else
                
                l_B1 = (l_CurrentWeigthedTristimuli' * l_CurrentWeigthedTristimuli)^(-1) * ...
                    l_CurrentWeigthedTristimuli' * l_ReferenceWeigthedTristimuli;
            
            end
                                   
            o_B = l_B1;
            o_resnormBEstimation = [];
            
            %% test other impl

            % base matlab lin problem: recommended, ls if over determined
            l_B2 = mldivide(i_CurrentTristimuli, i_ReferenceTristimuli);

            % base matlab lin problem: alternative, ls if over determined
            l_B3 = linsolve(i_CurrentTristimuli, i_ReferenceTristimuli);
            
            %lsqlin with constraints must be used separately
            l_B3a1 = lsqlin(i_CurrentTristimuli, i_ReferenceTristimuli(: , 1), ...
                [], [], [], [], [],[]);
            l_B3a2 = lsqlin(i_CurrentTristimuli, i_ReferenceTristimuli(: , 2), ...
                [], [], [], [], [],[]);
            l_B3a3 = lsqlin(i_CurrentTristimuli, i_ReferenceTristimuli(: , 3), ...
                [], [], [], [], [],[]);
            
            l_B3a = [l_B3a1, l_B3a2, l_B3a3];
            
            % optim toolbox: can be used also, although it is a linear problem
            if isempty(i_weightingVec)
                
                [l_B4, o_resnormBEstimation] = ...
                    lsqnonlin(@(l_B0) OFT_IDT_MeritFunctionCoreXYZ(l_B0, i_ReferenceTristimuli, i_CurrentTristimuli), l_B0);
            
            else
                
                [l_B4W, o_resnormBEstimationW] = ...
                    lsqnonlin(@(l_B0) OFT_IDT_MeritFunctionCoreXYZ_Weighted ...
                    (l_B0, l_ReferenceWeigthedTristimuli, i_CurrentTristimuli, l_W), l_B0);
            
            end
            
        otherwise
    end
    
    o_B(abs(o_B) < 0.0001) = 0.0;

end

%%CAM Models

function o_F = OFT_IDT_MeritFunctionCoreLuv(i_B0, i_ReferenceTristimuli, i_CurrentTristimuli, i_M, i_w)

    l_B = zeros(3, 3);
    i_B = i_B0;
    
    l_B(1, 1) = i_B(1, 1);
    l_B(2, 1) = i_B(2, 1);
    l_B(3, 1) = i_B(3, 1);
    
    l_B(1, 2) = i_B(1, 2);
    l_B(2, 2) = i_B(2, 2);
    l_B(3, 2) = i_B(3, 2);
    
    l_B(1, 3) = 1 - i_B(1, 1) - i_B(1, 2);
    l_B(2, 3) = 1 - i_B(2, 1) - i_B(2, 2);
    l_B(3, 3) = 1 - i_B(3, 1) - i_B(3, 2);

    k = 1 : size(i_ReferenceTristimuli,2);

    o_F = i_ReferenceTristimuli(:,k) -...
            HDM_OFT_ColorConversions.OFT_CIELuv(i_M * l_B * i_CurrentTristimuli(:,k), i_w);

    l_norm = sqrt(sum(abs(o_F).^2, 1));

    o_F = l_norm';

end

function o_F = OFT_IDT_MeritFunctionCoreLuv_Weighted(i_B0, i_ReferenceTristimuli, i_CurrentTristimuli, i_M, i_w, i_W)

    l_B = zeros(3, 3);
    i_B = i_B0;
    
    l_B(1, 1) = i_B(1, 1);
    l_B(2, 1) = i_B(2, 1);
    l_B(3, 1) = i_B(3, 1);
    
    l_B(1, 2) = i_B(1, 2);
    l_B(2, 2) = i_B(2, 2);
    l_B(3, 2) = i_B(3, 2);
    
    l_B(1, 3) = 1 - i_B(1, 1) - i_B(1, 2);
    l_B(2, 3) = 1 - i_B(2, 1) - i_B(2, 2);
    l_B(3, 3) = 1 - i_B(3, 1) - i_B(3, 2);
    
    k = 1 : size(i_ReferenceTristimuli,2);

    o_F = i_ReferenceTristimuli(:,k) -...
            (i_W * HDM_OFT_ColorConversions.OFT_CIELuv(i_M * l_B * i_CurrentTristimuli(:,k), i_w)')';

    l_norm = sqrt(sum(abs(o_F).^2, 1));

    o_F = l_norm';

end

function o_F = OFT_IDT_MeritFunctionCoreLab(i_B0, i_ReferenceTristimuli, i_CurrentTristimuli, i_M, i_w)

    l_B = zeros(3, 3);
    i_B = i_B0;
    
    l_B(1, 1) = i_B(1, 1);
    l_B(2, 1) = i_B(2, 1);
    l_B(3, 1) = i_B(3, 1);
    
    l_B(1, 2) = i_B(1, 2);
    l_B(2, 2) = i_B(2, 2);
    l_B(3, 2) = i_B(3, 2);
    
    l_B(1, 3) = 1 - i_B(1, 1) - i_B(1, 2);
    l_B(2, 3) = 1 - i_B(2, 1) - i_B(2, 2);
    l_B(3, 3) = 1 - i_B(3, 1) - i_B(3, 2);

    k = 1 : size(i_ReferenceTristimuli,2);

    o_F = i_ReferenceTristimuli(:,k) - ...
            HDM_OFT_ColorConversions.OFT_CIELab(i_M * l_B * i_CurrentTristimuli(:,k)', i_w)';

    % p2 norm per patch, it is also delta E;  ?or delta E 2000 
    o_F = sqrt(sum(o_F.^2, 2)); 
    
    % additional sqrt, cause client will square it, so overall it is ACES
    % method
    l_norm = sqrt(o_F);

    o_F = l_norm;
    
    %% you can test delta E 2000, but convergence seems to be bad
    
%     l_dE2000 = sqrt(deltaE2000(i_ReferenceTristimuli(:,k), ...
%         HDM_OFT_ColorConversions.OFT_CIELab(i_M * l_B * i_CurrentTristimuli(:,k)', i_w)'));
%     
%     % additional sqrt, cause client will square it, so overall it is ACES
%     % method
%     l_norm = sqrt(l_dE2000');
% 
%     o_F = l_norm;

end

function o_F = OFT_IDT_MeritFunctionCoreLab_Weighted(i_B0, i_ReferenceTristimuli, i_CurrentTristimuli, i_M, i_w, i_W)

    l_B = zeros(3, 3);
    i_B = i_B0;
    
    l_B(1, 1) = i_B(1, 1);
    l_B(2, 1) = i_B(2, 1);
    l_B(3, 1) = i_B(3, 1);
    
    l_B(1, 2) = i_B(1, 2);
    l_B(2, 2) = i_B(2, 2);
    l_B(3, 2) = i_B(3, 2);
    
    l_B(1, 3) = 1 - i_B(1, 1) - i_B(1, 2);
    l_B(2, 3) = 1 - i_B(2, 1) - i_B(2, 2);
    l_B(3, 3) = 1 - i_B(3, 1) - i_B(3, 2);

    k = 1 : size(i_ReferenceTristimuli,2);

    o_F = i_ReferenceTristimuli(:,k) - ...
            i_W * HDM_OFT_ColorConversions.OFT_CIELab(i_M * l_B * i_CurrentTristimuli(:,k)', i_w)';

    l_norm = sqrt(sum(abs(o_F).^2, 1)); % or delta E 2000

    o_F = l_norm';

end

%% Primary Models

function o_F = OFT_IDT_MeritFunctionCoreXYZ(i_B0, i_ReferenceTristimuli, i_CurrentTristimuli)

    l_B = zeros(3, 3);
    i_B = i_B0;
    
    l_B(1, 1) = i_B(1, 1);
    l_B(2, 1) = i_B(2, 1);
    l_B(3, 1) = i_B(3, 1);
    
    l_B(1, 2) = i_B(1, 2);
    l_B(2, 2) = i_B(2, 2);
    l_B(3, 2) = i_B(3, 2);
    
    l_B(1, 3) = 1 - i_B(1, 1) - i_B(1, 2);
    l_B(2, 3) = 1 - i_B(2, 1) - i_B(2, 2);
    l_B(3, 3) = 1 - i_B(3, 1) - i_B(3, 2);

    o_F = i_ReferenceTristimuli - i_CurrentTristimuli * l_B;

end

function o_F = OFT_IDT_MeritFunctionCoreXYZ_Weighted(i_B0, i_ReferenceTristimuliWeighted, i_CurrentTristimuli, i_W)

    l_B = zeros(3, 3);
    i_B = i_B0;
    
    l_B(1, 1) = i_B(1, 1);
    l_B(2, 1) = i_B(2, 1);
    l_B(3, 1) = i_B(3, 1);
    
    l_B(1, 2) = i_B(1, 2);
    l_B(2, 2) = i_B(2, 2);
    l_B(3, 2) = i_B(3, 2);
    
    l_B(1, 3) = 1 - i_B(1, 1) - i_B(1, 2);
    l_B(2, 3) = 1 - i_B(2, 1) - i_B(2, 2);
    l_B(3, 3) = 1 - i_B(3, 1) - i_B(3, 2);

    o_F = i_ReferenceTristimuliWeighted - i_W * i_CurrentTristimuli * l_B;

end

