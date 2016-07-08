function out = HDM_OFT_GetSpectralResponse(OFT_SpectralResponse)

    OFT_Env=HDM_OFT_InitEnvironment();

    if(exist('OFT_SpectralResponse','var') == 0)   
        
        disp('using default spectral response');
        OFT_SpectralResponse = 'RICD'
        
    else   
        
        default = 0;   
        disp('using spectral response');
        OFT_SpectralResponse
        
    end

    %% read
    
    if(isempty(strfind(OFT_SpectralResponse,'.')))
        
        switch OFT_SpectralResponse
            case 'RICD'
                
                OFT_RICDSpectralResponseFile = strcat(OFT_Env.OFT_RootDataDir,'/cameraResponsesReference/ricd.xlsx');
                [ndata, text, alldata] =xlsread(OFT_RICDSpectralResponseFile);
                OFT_CameraSpectralResponse_1nm_CIE31Range = ndata;
                
            otherwise
        end
        
    else
        
     	[l_p,l_n,l_ext]=fileparts(OFT_SpectralResponse);
        if(strcmp(l_ext,'.xls') || strcmp(l_ext,'.xlsx'))
            
            [ndata, text, alldata] =xlsread(OFT_SpectralResponse);
            OFT_CameraSpectralResponse_1nm_CIE31Range = ndata;
            
        elseif(strcmp(l_ext,'.csv'))
        
            OFT_CameraSpectralResponse=csvread(OFT_SpectralResponse);
            OFT_CameraSpectralResponse_Wavelength=OFT_CameraSpectralResponse(:,1)';

            OFT_CameraSpectralResponse_Rrelative=OFT_CameraSpectralResponse(:,2)';

            if(size(OFT_CameraSpectralResponse,2)==4)
                OFT_CameraSpectralResponse_Grelative=OFT_CameraSpectralResponse(:,3)';
                OFT_CameraSpectralResponse_Brelative=OFT_CameraSpectralResponse(:,4)';
            else
                OFT_CameraSpectralResponse_G1relative=OFT_CameraSpectralResponse(:,3)';
                OFT_CameraSpectralResponse_G2relative=OFT_CameraSpectralResponse(:,4)';
                OFT_CameraSpectralResponse_Grelative=0.5 * (OFT_CameraSpectralResponse_G1relative+OFT_CameraSpectralResponse_G2relative);
                OFT_CameraSpectralResponse_Brelative=OFT_CameraSpectralResponse(:,5)';
            end

            OFT_CameraSpectralResponse=...
                [OFT_CameraSpectralResponse_Wavelength;
                OFT_CameraSpectralResponse_Rrelative;
                OFT_CameraSpectralResponse_Grelative;
                OFT_CameraSpectralResponse_Brelative];

            OFT_CameraSpectralResponse_1nm_CIE31Range=...
                [360:1:830;
                interp1(OFT_CameraSpectralResponse(1,:),OFT_CameraSpectralResponse(2,:),360:1:830,'pchip',0);
                interp1(OFT_CameraSpectralResponse(1,:),OFT_CameraSpectralResponse(3,:),360:1:830,'pchip',0);
                interp1(OFT_CameraSpectralResponse(1,:),OFT_CameraSpectralResponse(4,:),360:1:830,'pchip',0)];
        
        end
        
    end
    
    %% normalize to green maximum
    
    l_greenMax = max(OFT_CameraSpectralResponse_1nm_CIE31Range(3,:));
    
    OFT_CameraSpectralResponse_1nm_CIE31Range(2, :) = OFT_CameraSpectralResponse_1nm_CIE31Range(2, :) ./ l_greenMax;
    OFT_CameraSpectralResponse_1nm_CIE31Range(3, :) = OFT_CameraSpectralResponse_1nm_CIE31Range(3, :) ./ l_greenMax;
    OFT_CameraSpectralResponse_1nm_CIE31Range(4, :) = OFT_CameraSpectralResponse_1nm_CIE31Range(4, :) ./ l_greenMax;

    %% show
    if(exist('default','var')==0)
        
        figure
        plot(OFT_CameraSpectralResponse_1nm_CIE31Range(1,:),OFT_CameraSpectralResponse_1nm_CIE31Range(4,:),...
            OFT_CameraSpectralResponse_1nm_CIE31Range(1,:),OFT_CameraSpectralResponse_1nm_CIE31Range(3,:),...
            OFT_CameraSpectralResponse_1nm_CIE31Range(1,:),OFT_CameraSpectralResponse_1nm_CIE31Range(2,:));
        xlabel('wavelength in nm')
        ylabel('green maximum normalized relative spectral response')
        legend({'b','g','r'})
        grid on
        grid minor
        title('spectral response of camera system');
    
    end

    out=OFT_CameraSpectralResponse_1nm_CIE31Range;

end

