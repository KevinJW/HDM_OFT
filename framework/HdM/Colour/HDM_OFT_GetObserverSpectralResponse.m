function o_observerSpectralResponse = HDM_OFT_GetObserverSpectralResponse(i_observerSpectralResponse)

    l_res = strfind(i_observerSpectralResponse,'.');

    if(isempty(l_res))

        o_observerSpectralResponse = HDM_OFT_CIEStandard.GetStandardObserverCurves(i_observerSpectralResponse);

    else

        l_rawResponse = csvread(i_observerSpectralResponse)';
        
        l_nmBase = 360:1:830;
        l_responseR = interp1(l_rawResponse(1, :), l_rawResponse(2, :), l_nmBase,'pchip',0);
        l_responseG = interp1(l_rawResponse(1, :), l_rawResponse(3, :), l_nmBase,'pchip',0);
        l_responseB = interp1(l_rawResponse(1, :), l_rawResponse(4, :), l_nmBase,'pchip',0);
        
        o_observerSpectralResponse = [l_nmBase; l_responseR; l_responseG; l_responseB];

    end

end

