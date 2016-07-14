classdef HDM_OFT_SpectrumExportImport
    methods(Static)
        function o_spectrum=ImportSpectrum(i_spectrumFile)

            %% defaults
            l_Env=HDM_OFT_InitEnvironment(); 

            HDM_OFT_Utils.OFT_DispTitle('import spectrum');
            
            %% uprtek
            
            try 
                
            l_fid = fopen(i_spectrumFile);
            l_out = textscan(l_fid,'%s%s','delimiter','\t');
            fclose(l_fid);

            l_s = size(l_out);
            
            if(l_s(1) == 1 && l_s(2) == 2)
                            
                l_c1 = l_out(1, 1);
                l_c1 = l_c1{1, 1};
                l_c2 = l_out(1, 2);
                l_c2 = l_c2{1, 1};
                
                l_w = [];
                l_i = [];
                
                for cur = 1 : size(l_c1, 1)
                    
                    l_str = l_c1 {cur, 1};
                    
                    if(findstr(l_str,'nm'))
                        
                        l_w = [l_w, str2num(strrep(l_str, 'nm', ''))];
                        l_i = [l_i, str2double(l_c2(cur))];
                        
                    end
                    
                end  
                
                o_spectrum = [l_w; l_i];
                
                l_nmBase=360:1:830;
                o_spectrum = [l_nmBase; interp1(o_spectrum(1,:), o_spectrum(2,:),l_nmBase,'pchip',0)];
                
                return;
                            
            end
            
            catch
            
            ;
                
            end
            
            [l_p,l_n,l_ext]=fileparts(i_spectrumFile);
            
            %% Open Film Tools spectral database
            if(strcmp(l_ext,'.csv'))
                try

                    l_table = HDM_OFT_SpectrumExportImport.read_mixed_csv(i_spectrumFile, ';');

                    l_table(1, :) = cellfun(@(s) {strrep(s, 'nm', '')},l_table(1, :));
                    l_table(1, :) = cellfun(@(s) {str2num(s)},l_table(1, :));
                    l_table(2, :) = cellfun(@(s) {str2double(s)},l_table(2, :));

                    o_spectrum = cell2mat([l_table(1, 2 : size(l_table, 2)); l_table(2, 2 : size(l_table, 2))]);

                    if(isnan(o_spectrum(2,1)))%%give third line a chance
                        
                        l_table(3, :) = cellfun(@(s) {str2double(s)},l_table(3, :));
                        o_spectrum = cell2mat([l_table(1, 2 : size(l_table, 2)); l_table(3, 2 : size(l_table, 2))]);
                        
                    end
                    
                    l_nmBase=360:1:830;
                    o_spectrum = [l_nmBase; interp1(o_spectrum(1,:), o_spectrum(2,:),l_nmBase,'pchip',0)];
                    return;

                catch

                    ;

                end
            end
                      
            %% EyeOne
                        
            if(strcmp(l_ext,'.xls') || strcmp(l_ext,'.xlsx'))
                
                try
                    
                [ndata, text, alldata] =xlsread(i_spectrumFile);
                
                if isempty(text)
                    
                    o_spectrum=ndata(1:2,:);
                    
                else
                    
                    waveLength=strrep(text(1,:), 'nm', '');
                    wvSize=size(waveLength);
                    waveLength2=str2double(waveLength(2:wvSize(2)));
                    o_spectrum = [waveLength2;ndata];
                    
                end
                
                l_nmBase=360:1:830;
                o_spectrum = [l_nmBase; interp1(o_spectrum(1,:), o_spectrum(2,:),l_nmBase,'pchip',0)];
                
                catch %%csv mimics excel but csv with non numeric data, which matlab does not support via csvread
                    
                    l_table = HDM_OFT_SpectrumExportImport.read_mixed_csv(i_spectrumFile, '\t');
                    
                    l_table(1, :) = cellfun(@(s) {strrep(s, 'nm', '')},l_table(1, :));
                    l_table(1, :) = cellfun(@(s) {str2num(s)},l_table(1, :));
                    l_table(2, :) = cellfun(@(s) {strrep(s, '.', '')},l_table(2, :));%//!!!
                    l_table(2, :) = cellfun(@(s) {str2double(s)},l_table(2, :));
                    
                    o_spectrum = cell2mat([l_table(1, 2 : size(l_table, 2)); l_table(2, 2 : size(l_table, 2))]);
                    
                    l_nmBase=360:1:830;
                    o_spectrum = [l_nmBase; interp1(o_spectrum(1,:), o_spectrum(2,:),l_nmBase,'pchip',0)];
                    
                end
                
            else        
                
                l_csvData = csvread(i_spectrumFile);
                                                              
            end
            
            
        end
        
        %%http://stackoverflow.com/questions/4747834/import-csv-file-with-mixed-data-types
        function lineArray = read_mixed_csv(fileName,delimiter)
          fid = fopen(fileName,'r');   %# Open the file
          lineArray = cell(100,1);     %# Preallocate a cell array (ideally slightly
                                       %#   larger than is needed)
          lineIndex = 1;               %# Index of cell to place the next line in
          nextLine = fgetl(fid);       %# Read the first line from the file
          while ~isequal(nextLine,-1)         %# Loop while not at the end of the file
            lineArray{lineIndex} = nextLine;  %# Add the line to the cell array
            lineIndex = lineIndex+1;          %# Increment the line index
            nextLine = fgetl(fid);            %# Read the next line from the file
          end
          fclose(fid);                 %# Close the file
          lineArray = lineArray(1:lineIndex-1);  %# Remove empty cells, if needed
          for iLine = 1:lineIndex-1              %# Loop over lines
            lineData = textscan(lineArray{iLine},'%s',...  %# Read strings
                                'Delimiter',delimiter);
            lineData = lineData{1};              %# Remove cell encapsulation
            if strcmp(lineArray{iLine}(end),delimiter)  %# Account for when the line
              lineData{end+1} = '';                     %#   ends with a delimiter
            end
            lineArray(iLine,1:numel(lineData)) = lineData;  %# Overwrite line data
          end
        end
        
    end
end