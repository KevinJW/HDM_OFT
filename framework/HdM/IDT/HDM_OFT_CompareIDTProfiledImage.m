function HDM_OFT_CompareIDTProfiledImage(OFT_In_ReferenceImage2View, OFT_In_TransformedImage2View)

OFT_Env=HDM_OFT_InitEnvironment();

HDM_OFT_Utils.OFT_DispTitle('compare transformed image');

outputVideo = VideoWriter(fullfile(OFT_Env.OFT_ProcessPath,'Eval2Reference.avi'), 'Uncompressed AVI');
outputVideo.FrameRate = 24;
open(outputVideo);

OFT_TransformedImage2View = imresize(OFT_In_TransformedImage2View, [NaN 1000]);
OFT_ReferenceImage2View = imresize(OFT_In_ReferenceImage2View, [size(OFT_TransformedImage2View , 1) 1000]);

oft_pixNorm=1.0/max(OFT_TransformedImage2View(:));
oft_pixNormRef=1.0/max(OFT_ReferenceImage2View(:));

OFT_ImageOriginalDouble1Base=oft_pixNorm*double(OFT_TransformedImage2View);
OFT_ImageOriginalOrgDouble1Base=oft_pixNormRef*double(OFT_ReferenceImage2View);

for ii = 1:2
    for pulse=1.0:-0.1:0
        l_sum = pulse*OFT_ImageOriginalOrgDouble1Base+(1-pulse)*OFT_ImageOriginalDouble1Base;    
        l_sum(l_sum > 1.0) = 1.0;
        l_sum(l_sum < 0.0) = 0.0;
        writeVideo(outputVideo, l_sum);
    end
    for pulse=0:0.1:1.0
        l_sum = pulse*OFT_ImageOriginalOrgDouble1Base+(1-pulse)*OFT_ImageOriginalDouble1Base;
        l_sum(l_sum > 1.0) = 1.0;
        l_sum(l_sum < 0.0) = 0.0;
        writeVideo(outputVideo, l_sum);
    end
end

close(outputVideo);

implay(fullfile(OFT_Env.OFT_ProcessPath,'Eval2Reference.avi'));

HDM_OFT_Utils.OFT_DispTitle('compare transformed image succesfully finished');

end
