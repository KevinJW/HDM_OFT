function o_Env_UI=HDM_OFT_InitEnvironment_UI...
    (i_DoPlot, ...
    i_OutPlotColourDir, ...
    i_OutPlotBWDir, ...
    i_SavePlot, ...
    i_ExportPlot, ...
    i_BW)

    global globalOFT_Env_UI;

    if(isempty(globalOFT_Env_UI))
       
    close all;

    % The properties we've been using in the figures

    % Defaults
    
    %                 Default     Paper     Presentation
    % Width           5.6         varies	varies
    % Height          4.2         varies	varies
    % AxesLineWidth   0.5         0.75      1
    % FontSize        10          8         14
    % LineWidth       0.5         1.5       2
    % MarkerSize      6           8         12
      
    l_width = 4;     % Width in inches
    l_height = 3;    % Height in inches
    l_alw = 0.75;    % AxesLineWidth
    l_fsz = 10;      % Fontsize
    l_lw = 1.5;      % LineWidth
    l_msz = 8;       % MarkerSize
    
    set(0,'DefaultLineLineWidth', l_lw);   % set the default line width to lw
    set(0,'DefaultLineMarkerSize', l_msz); % set the default line marker size to msz
    set(0,'DefaultLineLineWidth', l_lw);   % set the default line width to lw
    set(0,'DefaultLineMarkerSize', l_msz); % set the default line marker size to msz

    % Set the default Size for display
    defpos = get(0,'DefaultFigurePosition');
    
    l_screenSize = get(0,'screensize');
    set(0,'DefaultFigurePosition', [30 (l_screenSize(4) - 30) l_width * 100, l_height * 100]);
    
    % Set the defaults for saving/printing to a file
    set(0,'DefaultFigureInvertHardcopy','on'); % This is the default anyway
    set(0,'DefaultFigurePaperUnits','inches'); % This is the default anyway
    l_defsize = get(gcf, 'PaperSize');
    l_left = (l_defsize(1)- l_width)/2;
    l_bottom = (l_defsize(2)- l_height)/2;
    l_defsize = [l_left, l_bottom, l_width, l_height];
    set(0, 'DefaultFigurePaperPosition', l_defsize);

    % Change default axes fonts.
    set(0,'DefaultAxesFontName', 'Arial');
    set(0,'DefaultAxesFontSize', l_fsz);

    % Change default text fonts.
    set(0,'DefaultTextFontname', 'Arial');
    set(0,'DefaultTextFontSize', l_fsz);
    
    % default grid on
    set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on')

    %set(0,'DefaultAxesLineStyleOrder','-|--|:|-.|.');

    % dco=get(0,'DefaultAxesColorOrder');
    
    l_mainVal = 0.7;
    l_sideVal = 0.0;
    
    set(0,...
        'DefaultAxesLineStyleOrder','-|--|-.|-..|.|--.|:',...
        'DefaultAxesColorOrder',...
    ...%    [0.0 0.0 0.0;0.45 0.45 0.45;0.7 0.7 0.7]); % grayish
            [l_mainVal, l_sideVal, l_sideVal;... %rgb cmy ordered
            l_sideVal, l_mainVal, l_sideVal;...
            l_sideVal, l_sideVal, l_mainVal]);
%             l_mainVal, l_sideVal, l_mainVal;...
%             l_mainVal, l_mainVal, l_sideVal;...
%             l_sideVal, l_mainVal, l_mainVal]);    

    l_Env_UI = HDM_OFT_Environment_UI();
    
    l_Env_UI.OFT_DoPlot = i_DoPlot;

    l_Env_UI.OFT_PlotColourDir = i_OutPlotColourDir;
    l_Env_UI.OFT_PlotBWDir = i_OutPlotBWDir;
    
    l_Env_UI.OFT_SavePlot = i_SavePlot;
    l_Env_UI.OFT_ExportPlotToClipBoard = i_ExportPlot;
    
    l_Env_UI.OFT_BW = i_BW;

    o_Env_UI=l_Env_UI;
    globalOFT_Env_UI=l_Env_UI;

    %clear classes %commented out due to in arg preserving
    delete(findall(0,'Type','figure'));

    else
        o_Env_UI=globalOFT_Env_UI;    
    end

end
        

        
