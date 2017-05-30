function HDM_OFT_UI_ImShowAndSave...
    (i_data, i_title)

    l_lastPos = get(0,'DefaultFigurePosition');
        
    l_h =  findobj('type','figure');
    l_n = length(l_h);

    if l_n > 0
        
        l_lastFig = gcf;
    
        l_n = get(gcf,'Number');  
    
        l_lastPos = get(l_lastFig,'Position');
        
    end   
    
    l_defpos = get(0,'DefaultFigurePosition');
    
    l_shift = 50;
    
    l_screenSize = get(0,'screensize');
    
    l_dy = 800 - l_n * l_shift;

    l_size = size(i_data, 1);

    l_figure = figure('Name', i_title, 'Position', [l_lastPos(1) + l_shift, l_dy, l_defpos(3), l_defpos(4)]); clf;
    
    if iscell(i_data)%{1,1}(1, 3)
        
        for cur = 1 : l_size

                l_t = i_data{cur,1}(1, 1);
                l_t2 = i_data{cur,1}(1, 2);
                l_s = i_data{cur,1}(1, 3);
                
                plot(l_t{1, 1}, l_t2{1, 1}, l_s{1, 1});

            hold on

        end
    
    else
               
        for cur = 2 : l_size

            plot(i_data(1, :), i_data(cur, :));

            hold on

        end
        
        xlim([min(i_data(1, :)) max(i_data(1, :))]);
        
    end
    
    if(exist('i_legend','var'))
        
        legend(i_legend);
    
    end
    
    xlabel(i_xLabel);
    ylabel(i_yLabel);
    
    % ylim([-0.1 1.1])
    
    if(exist('i_textArgs','var'))
        
        for cur = 1 : size(i_textArgs, 1)
            
            l_curTextArg1 = i_textArgs{cur,1}(1, 1);
            l_curTextArg2 = i_textArgs{cur,1}(1, 2);
            l_curTextArg3 = i_textArgs{cur,1}(1, 3);
            l_curTextArg4 = i_textArgs{cur,1}(1, 4);
            l_curTextArg5 = i_textArgs{cur,1}(1, 5);
            l_curTextArg6 = i_textArgs{cur,1}(1, 6);
            l_curTextArg7 = i_textArgs{cur,1}(1, 7);
            
            text(l_curTextArg1{1, 1}, l_curTextArg2{1, 1}, ...
                l_curTextArg3{1, 1}, ...
                l_curTextArg4{1, 1}, l_curTextArg5{1, 1}, l_curTextArg6{1, 1}, l_curTextArg7{1, 1})
            
        end
              
    end

    l_name = get(l_figure, 'name');
    
    l_shortFileName = strcat('Figure', {' '}, num2str(l_n), {' '}, l_name);
    l_shortFileName = l_shortFileName{1, 1};
    
    global globalOFT_Env_UI;
    if (globalOFT_Env_UI.OFT_ExportPlotToClipBoard)  
        print(l_figure, '-clipboard', '-dbitmap');
    end

    if (globalOFT_Env_UI.OFT_SavePlot)
        if(globalOFT_Env_UI.OFT_BW)
            savefig(l_figure, strcat(globalOFT_Env_UI.OFT_PlotBWDir, '/', l_shortFileName),'compact')
            print(l_figure, strcat(globalOFT_Env_UI.OFT_PlotBWDir, '/', l_shortFileName),'-dpng','-r300');
        else
            savefig(l_figure, strcat(globalOFT_Env_UI.OFT_PlotColourDir, '/', l_shortFileName),'compact');
            print(l_figure, strcat(globalOFT_Env_UI.OFT_PlotColourDir, '/', l_shortFileName),'-dpng','-r300');
        end
    end
    
end