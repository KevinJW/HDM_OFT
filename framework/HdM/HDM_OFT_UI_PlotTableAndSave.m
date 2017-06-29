function HDM_OFT_UI_PlotTableAndSave...
    (i_data, i_title, i_cNames, i_rNames)

    l_lastPos = get(0,'DefaultFigurePosition');
        
    l_h =  findobj('type','figure');
    l_n = length(l_h) + 1;

    if l_n > 1
        
        l_lastFig = gcf;
    
        l_n = get(gcf,'Number');  
    
        l_lastPos = get(l_lastFig,'Position');
        
    end   
    
    l_defpos = get(0,'DefaultFigurePosition');
    
    l_shift = 50;
    
    l_screenSize = get(0,'screensize');
    
    l_dy = 800 - l_n * l_shift;

    l_figure = figure('Name', i_title, 'Position', [l_lastPos(1) + l_shift, l_dy, l_defpos(3), l_defpos(4)]);
    l_t = uitable(l_figure, 'Data', i_data,...
                'ColumnName', i_cNames, 'RowName', i_rNames);%, ...
                %'Position', [0 0 100 100]);

    l_t.Position(3) = l_t.Extent(3);
    %l_t.Position(4) = l_t.Extent(4);
    
    l_name = get(l_figure, 'name');
    
    l_shortFileName = strcat('Figure', {' '}, num2str(l_n), {' '}, l_name, 'csv');
    l_shortFileName = l_shortFileName{1, 1};
    
    global globalOFT_Env_UI;
    if (globalOFT_Env_UI.OFT_ExportPlotToClipBoard)  
        
        clipboard('copy', i_data);
        
    end

    if (globalOFT_Env_UI.OFT_SavePlot)
        
        HDM_OFT_CSVWrite(l_shortFileName, i_data, i_cNames);
        
    end   
    
end