function o_IDTProfileFileName = HDM_OFT_WriteIDTProfileAndStatEntry...
    (i_MeasuredCameraResponseFileName, i_resnormBEstimation, i_IDT_B, i_IDT_b,...
    i_NeutralsCompensation, i_IlluminantSpectrum, i_IlluminantWhitePoint, i_ReferenceDomain, i_ErrorMinimizationDomain, ...
    i_meanDeltaE2000, i_stdDevDeltaE2000)

    l_Env = HDM_OFT_InitEnvironment();

    if size(i_IlluminantSpectrum, 1) == 2
        
        IlluminantStr=strcat('Illuminant_', 'fromData');
        
    elseif(isempty(strfind(i_IlluminantSpectrum,'.')))
        
        IlluminantStr=strcat('Illuminant_',i_IlluminantSpectrum);
        
    else
        
        [OFT_IlluminantSpectrumPath,OFT_IlluminantSpectrumName,OFT_IlluminantSpectrumExt] = fileparts(i_IlluminantSpectrum);
        IlluminantStr=strcat('Illuminant_FromFile_',OFT_IlluminantSpectrumName);
        
    end

    fin = fopen(strcat(l_Env.OFT_ConstraintsPath,'/IDT_Template.txt'));
    idtCreationDateStr=datestr(now,'yyyy-mm-dd_HH.MM.SS.FFF');
    OFT_IDT_File=strcat(l_Env.OFT_ProcessPath,'/IDT_',IlluminantStr,'_',idtCreationDateStr,'.ctl');
    fout = fopen(OFT_IDT_File,'wt');

    if ~exist(strcat(l_Env.OFT_StatisticsPath,'/IDTStat.csv'),'file')
        foutStat = fopen(strcat(l_Env.OFT_StatisticsPath,'/IDTStat.csv'),'wt');
        fprintf(foutStat,'idt file\t\t , measurement file\t\t , resnorm , B11 , B12 , B13 , B21 , B22 , B23 , B31 , B32 , B33 , scene adopted white, neutrals compensation, colour domain for error minimization, mean Delta E 2000, std dev Delta E 2000\n');
    else
        foutStat = fopen(strcat(l_Env.OFT_StatisticsPath,'/IDTStat.csv'),'at');
    end

    fprintf(foutStat,'%s\t , ', strcat('IDT_',IlluminantStr,'_',idtCreationDateStr,'.ctl'));
    
    if size(i_MeasuredCameraResponseFileName, 1) == 4
        
        fprintf(foutStat,'%s\t , ', 'from Data');
     
    else
        
        [~,OFT_MeasuredCameraResponseName,OFT_MeasuredCameraResponseExt] = fileparts(i_MeasuredCameraResponseFileName);
        fprintf(foutStat,'%s\t , ', strcat(OFT_MeasuredCameraResponseName,OFT_MeasuredCameraResponseExt));
    
    end
    
    fprintf(foutStat,'%e , ', i_resnormBEstimation);

    while ~feof(fin)
       S = fgetl(fin);
       %s = strrep(s, '118520', '118521');

       if(strfind(S, '%'))

            if(strfind(S, '%IDT_DATE%'))

               fprintf(fout,'// Creation Date: %s\n',idtCreationDateStr);

            elseif (strfind(S, '%IDT_ILLUMINANT%'))

               fprintf(fout,'// Illuminant %s\n', IlluminantStr);

            elseif(strfind(S, '%const float b['))

               fprintf(fout,'\tconst float b[] = { %f, %f, %f };\n',i_IDT_b(1),i_IDT_b(2),i_IDT_b(3));

            elseif(strfind(S, '%const float B1'))

               fprintf(fout,'\tconst float B[][] =     { { %f, %f, %f },\n',i_IDT_B(1,1),i_IDT_B(1,2),i_IDT_B(1,3));

               fprintf(foutStat,'%f , %f , %f , ',i_IDT_B(1,1),i_IDT_B(1,2),i_IDT_B(1,3));

            elseif(strfind(S, '%const float B2'))

               fprintf(fout,'\t\t\t  { %f, %f, %f },\n',i_IDT_B(2,1),i_IDT_B(2,2),i_IDT_B(2,3));

               fprintf(foutStat,'%f , %f , %f , ',i_IDT_B(2,1),i_IDT_B(2,2),i_IDT_B(2,3));

            elseif(strfind(S, '%const float B3'))

               fprintf(fout,'\t\t\t  { %f, %f, %f }};\n',i_IDT_B(3,1),i_IDT_B(3,2),i_IDT_B(3,3));

               fprintf(foutStat,'%f , %f , %f , ',i_IDT_B(3,1),i_IDT_B(3,2),i_IDT_B(3,3));

               fprintf(foutStat,'%s , %s , %s, %f, %f', IlluminantStr, i_NeutralsCompensation, i_ErrorMinimizationDomain, i_meanDeltaE2000, i_stdDevDeltaE2000);

           end

       else

            fprintf(fout,'%s\n',S);

       end
    end

    fprintf(foutStat,'\n');

    fclose(fin);
    fclose(fout);
    fclose(foutStat);
    
    %%write icc
    
    if 0 %HDM_OFT_IDT_ReferenceCamera.CIEType() == i_ReferenceDomain
        
        l_fileICC=strcat(l_Env.OFT_ProcessPath,'/IDT_', IlluminantStr,'_',idtCreationDateStr,'.icc ');
        l_bin = '/Volumes/Data/User/fuchs/Development/WriteICC/build-Write2ICC-Desktop_Qt_5_8_0_clang_64bit-Debug/Write2ICC ';

        l_IDT_b = i_IDT_b ./ min(i_IDT_b);
        l_IDT_b = l_IDT_b ./ l_IDT_b;

        %% scaling goes in lin curve max
        l_bStr = sprintf(' -b %f %f %f ',...
            l_IDT_b(1),l_IDT_b(2),l_IDT_b(3));

        [l_M_2xyz, l_wRef] = HDM_OFT_IDT_ReferenceCamera.GetDefinition(i_ReferenceDomain);
        
        l_XYZ4R = l_M_2xyz * i_IDT_B * ([1; 0; 0] .* (i_IDT_b ./ min(i_IDT_b)));
        l_XYZ4G = l_M_2xyz * i_IDT_B * ([0; 1; 0] .* (i_IDT_b ./ min(i_IDT_b)));
        l_XYZ4B = l_M_2xyz * i_IDT_B * ([0; 0; 1] .* (i_IDT_b ./ min(i_IDT_b)));
        
        l_M = [l_XYZ4R, l_XYZ4G, l_XYZ4B]';

        l_MStr = sprintf('-M %f %f %f %f %f %f %f %f %f',...
            l_M(1,1),l_M(1,2),l_M(1,3),...
            l_M(2,1),l_M(2,2),l_M(2,3),...
            l_M(3,1),l_M(3,2),l_M(3,3));

        l_w = i_IlluminantWhitePoint;

        l_w = l_w / l_w(2);

        l_whiteStr = sprintf(' -w %f %f %f ',...
            l_w(1),l_w(2),l_w(3));

        l_IlluminantWhitePoint = i_IlluminantWhitePoint / i_IlluminantWhitePoint(2);

        l_Mnc = HDM_OFT_ColorNeutralCompensations.OFT_GetNeutralCompensationMatrix(i_NeutralsCompensation, l_IlluminantWhitePoint, l_w);

        l_MncStr = sprintf('-Mnc %f %f %f %f %f %f %f %f %f',...
            l_Mnc(1,1),l_Mnc(1,2),l_Mnc(1,3),...
            l_Mnc(2,1),l_Mnc(2,2),l_Mnc(2,3),...
            l_Mnc(3,1),l_Mnc(3,2),l_Mnc(3,3));

        l_arg = [l_bin l_fileICC l_bStr l_MStr l_whiteStr l_MncStr];
        strcat(l_arg);

        system(l_arg);

        l_fileICC4OSX = strcat('~/Library/ColorSync/Profiles', '/IDT_', IlluminantStr,'_',idtCreationDateStr,'.icc');
        copyfile(l_fileICC, l_fileICC4OSX);


        %%

%         OFT_CIEStandardObserver_SpectralCurves=HDM_OFT_CIEStandard.GetStandardObserverCurves(HDM_OFT_CIEStandard.StandardObserver1931_2Degrees());
% 
%         OFT_Illuminant_Spectrum_1nm_CIE31Range=HDM_OFT_GetIlluminantSpectrum(i_IlluminantSpectrum);
% 
%         OFT_CameraSpectralResponse_1nm_CIE31Range = HDM_OFT_GetSpectralResponse(i_MeasuredCameraResponseFileName);
% 
%         [l_Tristimuli4Corner_TargetIlluminant,l_Chromaticities4Corner_TargetIlluminant]=...
%             HDM_OFT_TristimuliCreator.CreateFromSpectrum(...
%                     OFT_CIEStandardObserver_SpectralCurves,...
%                     OFT_Illuminant_Spectrum_1nm_CIE31Range,...
%                     OFT_CameraSpectralResponse_1nm_CIE31Range);  
% 
%         l_M_P = l_Tristimuli4Corner_TargetIlluminant ./ 100.0;
% 
%         l_M_P = l_M_P';
% 
%         l_M_PStr = sprintf('-M %f %f %f %f %f %f %f %f %f',...
%             l_M_P(1,1),l_M_P(1,2),l_M_P(1,3),...
%             l_M_P(2,1),l_M_P(2,2),l_M_P(2,3),...
%             l_M_P(3,1),l_M_P(3,2),l_M_P(3,3));
% 
%         l_fileICC_P=strcat(l_Env.OFT_ProcessPath,'/IDT_',IlluminantStr,'_',idtCreationDateStr,'_primaries.icc ');
% 
%         l_arg = [l_bin l_fileICC_P l_bStr l_M_PStr l_whiteStr l_MncStr];
%         strcat(l_arg);
% 
%         system(l_arg);
% 
%         l_fileICC4OSX = strcat('~/Library/ColorSync/Profiles', '/IDT_',IlluminantStr,'_',idtCreationDateStr,'_primaries.icc');
%         copyfile(l_fileICC_P, l_fileICC4OSX);

        %append to tiff example
    %     l_fid = fopen(l_fileICC);
    %     l_raw_profile_bytes = fread(l_fid, Inf, 'uint8=>uint8');
    %     fclose(l_fid);
    %     
    %     l_tif = Tiff('peppers.tif', 'r+');
    %     l_tif.setTag('ICCProfile', l_raw_profile_bytes);
    %     l_tif.rewriteDirectory();
    %     l_tif.close();

    end
    
    o_IDTProfileFileName=OFT_IDT_File;

end
