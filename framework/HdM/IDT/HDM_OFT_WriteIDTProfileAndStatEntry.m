function o_IDTProfileFileName = HDM_OFT_WriteIDTProfileAndStatEntry...
    (i_MeasuredCameraResponseFileName, i_resnormBEstimation, i_IDT_B, i_IDT_b,...
    i_NeutralsCompensation, i_IlluminantSpectrum, i_ErrorMinimizationDomain)

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
        fprintf(foutStat,'idt file\t\t , measurement file\t\t , resnorm , B11 , B12 , B13 , B21 , B22 , B23 , B31 , B32 , B33 , scene adopted white, neutrals compensation, colour domain for error minimization\n');
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

               fprintf(fout,'// Illuminant %s\n', i_IlluminantSpectrum);

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

               fprintf(foutStat,'%s , %s , %s', i_IlluminantSpectrum, i_NeutralsCompensation, i_ErrorMinimizationDomain);

           end

       else

            fprintf(fout,'%s\n',S);

       end
    end

    fprintf(foutStat,'\n');

    fclose(fin);
    fclose(fout);
    fclose(foutStat);
    
    o_IDTProfileFileName=OFT_IDT_File;

end
