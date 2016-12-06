function HDM_OFT_IDT_CreateBySpectralResponse_In_ZIPorXML(OFT_In_ClientData, OFT_In_ServerOutDir, OFT_In_TaskID)

    try

    %% init

    OFT_Env=HDM_OFT_InitEnvironment(OFT_In_TaskID);

    [OFT_AppProcessorPath,OFT_ProcessorName,OFT_ProcessorExt] = fileparts(mfilename('fullpath'));

    HDM_OFT_InitEnvironment_UI ...
        (true, strcat(OFT_AppProcessorPath,'/outPlotsInColour'), ...
        strcat(OFT_AppProcessorPath,'/outPlotsInBW'), ...
        false, false, false);    

    HDM_OFT_Utils.OFT_DispTitle('start task');
    HDM_OFT_Utils.OFT_DispSubTitle(strcat('ID: ',OFT_In_TaskID));
    HDM_OFT_Utils.OFT_DispSubTitle(strcat('client data: ',OFT_In_ClientData));
    HDM_OFT_Utils.OFT_DispSubTitle(strcat('output directory: ',OFT_In_ServerOutDir));

    OFT_LogFileName=strcat(OFT_In_ServerOutDir,'/',OFT_In_TaskID,'.status.xml');
    OFT_ProgressLogger=HDM_OFT_XML_Logger(OFT_LogFileName);
    OFT_ProgressLogger.LogUserMessage('start camera characterization');

    %% read client data

    [OFT_ClientDataPath,OFT_ClientDataName,OFT_ClientDataExt] = fileparts(OFT_In_ClientData);
    copyfile(OFT_ClientDataPath,OFT_Env.OFT_ProcessPath);
    IDTTaskData=HDM_OFT_IDT_PrepareClientData(strcat(OFT_Env.OFT_ProcessPath,'/',OFT_ClientDataName,OFT_ClientDataExt));
    IDTTaskData.ServerOutDir=OFT_In_ServerOutDir;

    OFT_ProgressLogger.LogUserMessage('estimate camera linearization');

    %% linearization, currently dummy, data at input must be linearized

    IDTTaskData.PreLinearisation_Out_LinCurve=HDM_OFT_CameraLinearization(IDTTaskData.PreLinearisation_In_LinFile);
    if(~strcmp(IDTTaskData.PreLinearisation_Out_LinCurve,''))
        copyfile(IDTTaskData.PreLinearisation_Out_LinCurve, OFT_In_ServerOutDir);   
    end

    OFT_ProgressLogger.LogUserMessage('estimate camera spectral response');


    % Decide use case: Spectral or color checker based method
    if isempty(IDTTaskData.SpectralResponse_In_LineCalibrationSpectrum) && isempty(IDTTaskData.SpectralResponse_In_LineCalibrationImage) && isempty(IDTTaskData.SpectralResponse_In_LightCalibrationSpectrum)
        IDTTaskData.SpectralResponse_Out_SpectralResponseFile=IDTTaskData.SpectralResponse_In_LightCalibrationImage;
    else

        %% get spectral response

        OFT_CameraResponse=HDM_OFT_CameraSpectralResponse(IDTTaskData);

        OFT_CameraReponseFile=strcat('/', OFT_In_TaskID, '_cameraResponse.csv');    
        IDTTaskData.SpectralResponse_Out_SpectralResponseFile=strcat(OFT_Env.OFT_ProcessPath,'/',OFT_CameraReponseFile);
        csvwrite(IDTTaskData.SpectralResponse_Out_SpectralResponseFile, OFT_CameraResponse');%//!!!why here
        copyfile(IDTTaskData.SpectralResponse_Out_SpectralResponseFile, OFT_In_ServerOutDir);   

        disp(strcat('camera spectral response file: ',IDTTaskData.SpectralResponse_Out_SpectralResponseFile));
    end

    %% generate profiles

    OFT_ProgressLogger.LogUserMessage('compute IDT profiles');    

    IDTTaskData.IDTCreationConstraints_Out_IDTFiles=HDM_OFT_IDT_ProfilesGeneration(IDTTaskData);

    for i=1:size(IDTTaskData.IDTCreationConstraints_Out_IDTFiles,2)
        copyfile(IDTTaskData.IDTCreationConstraints_Out_IDTFiles{i}, OFT_In_ServerOutDir);
    end

    OFT_ProgressLogger.LogUserMessage('finish camera characterization');

    HDM_OFT_Utils.OFT_DispTitle('task finished');
    HDM_OFT_Utils.OFT_DispSubTitle(strcat('ID: ',OFT_In_TaskID));
    HDM_OFT_Utils.OFT_DispSubTitle(strcat('client data: ',OFT_In_ClientData));
    HDM_OFT_Utils.OFT_DispSubTitle(strcat('output directory: ',OFT_In_ServerOutDir));
    commandwindow;

    catch err
        commandwindow;
        disp(getReport(err));
        OFT_ProgressLogger.LogUserMessage('error during camera characterization');
    end

end
        

        
