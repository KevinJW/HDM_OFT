classdef HDM_OFT_CIEStandard
    methods(Static)

        function [o_StandardObserverCurves] = ...
                GetStandardObserverCurves(OFT_In_StandardObserverType)
        
            
        if exist('OFT_In_StandardObserverType','var') && size(OFT_In_StandardObserverType, 1) == 4
    
            o_StandardObserverCurves = OFT_In_StandardObserverType;
            
            return;

        end    
            
        OFT_Env=HDM_OFT_InitEnvironment();    
        
        if(exist('OFT_In_StandardObserverType','var')==0)   
            disp('using default observer');
            OFT_StandardObserverType=HDM_OFT_CIEStandard.StandardObserver1931_2Degrees();
        else    
            default=0;
            disp('using given observer');
            OFT_StandardObserverType=OFT_In_StandardObserverType;
        end
            
        switch OFT_StandardObserverType

        case HDM_OFT_CIEStandard.StandardObserver1931_2Degrees()
            
            HDM_OFT_Utils.OFT_DispSubTitle('setup CIE 31 curves');
            %1 nm resolution
            OFT_CIE31_SpectralCurve=csvread(strcat(OFT_Env.OFT_ConstraintsPath,'/ciexyz31_1.csv'));
            OFT_CIE31_WaveLength=OFT_CIE31_SpectralCurve(:,1);
            OFT_CIE31_SpectralCurve_x=[OFT_CIE31_SpectralCurve(:,2)];
            OFT_CIE31_SpectralCurve_y=[OFT_CIE31_SpectralCurve(:,3)];
            OFT_CIE31_SpectralCurve_z=[OFT_CIE31_SpectralCurve(:,4)];
            
            OFT_WaveLengthRange=OFT_CIE31_WaveLength';
            OFT_StandardObserver_x=OFT_CIE31_SpectralCurve_x';
            OFT_StandardObserver_y=OFT_CIE31_SpectralCurve_y';
            OFT_StandardObserver_z=OFT_CIE31_SpectralCurve_z';
            
            o_StandardObserverCurves=[OFT_WaveLengthRange;...
                                            OFT_StandardObserver_x;OFT_StandardObserver_y;OFT_StandardObserver_z];

        otherwise

            HDM_OFTP_Utils.OFT_DispSubTitle('not implemented');

        end        
        
        end
 
        function out = StandardObserver1931_2Degrees()
            out='Standard Observer 1931 2 Degrees';
        end
      
        function [OFT_CIE31_colorValuePartsWaveLength,OFT_CIE31_colorValueParts_x,OFT_CIE31_colorValueParts_y,OFT_CIE31_colorValueParts_z]=...
                ColorValuePartsForSpectralColoursCurve(OFT_StandardObserverType)
           
            OFT_Env=HDM_OFT_InitEnvironment(); 
            
            if(~(OFT_StandardObserverType==HDM_OFT_CIEStandard.StandardObserver1931_2Degrees))
                return;
            end
            
            OFT_CIE31_colorValueParts=csvread(strcat(OFT_Env.OFT_ConstraintsPath,'/cccie31.csv'));
            OFT_CIE31_colorValuePartsWaveLength=OFT_CIE31_colorValueParts(:,1);
            OFT_CIE31_colorValueParts_x=[OFT_CIE31_colorValueParts(:,2)];
            OFT_CIE31_colorValueParts_y=[OFT_CIE31_colorValueParts(:,3)];
            OFT_CIE31_colorValueParts_z=[OFT_CIE31_colorValueParts(:,4)];

        end
    end
end