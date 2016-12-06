function HDM_OFT_CSVWrite...
    (i_shortFileName, i_data, i_header)

    global globalOFT_App_Env;
    
    l_fullFileName = strcat(globalOFT_App_Env.OFT_ProcessPath, '/', i_shortFileName);

    l_header = i_header; %dummy header
    l_commaHeader = [l_header; repmat({','}, 1, numel(l_header))]; %insert commaas
    l_commaHeader = l_commaHeader(:)';
    l_textHeader = cell2mat(l_commaHeader); %cHeader in text with commas

    %write header to file
    fid = fopen(l_fullFileName, 'w'); 
    fprintf(fid,'%s\n', l_textHeader);
    fclose(fid);

    %write data to end of file
    dlmwrite(l_fullFileName, i_data, '-append');
    
end