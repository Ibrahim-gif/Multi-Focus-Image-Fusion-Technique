% After all components done, combine code here

function outImg = main(id)
    % global var
    SML_N = 1;
    SML_T = 0;
    SML_s = 1;
    
    type = "color";
    if (id(1) == 'g')
        type = "grayscale";
    elseif (id(1) == 'n')
        type = "noisy";
    end

    % part 1, import and normalize
    original1 = im2double(imread(sprintf('testImg/%s/%s_1.tif', type, id)));
    original2 = im2double(imread(sprintf('testImg/%s/%s_2.tif', type, id)));

    % part 2, highpass lowpass split
    [LForiginal1, HForiginal1, LForiginal2, HForiginal2] = ...
        CutImage(original1, original2);
    % NOTE hf images are enhanced
    imwrite(LForiginal1, sprintf("outputImg/cut/%s_1_lf.tif", id));
    imwrite(HForiginal1 * 10, sprintf("outputImg/cut/%s_1_hf.tif", id));
    imwrite(LForiginal2, sprintf("outputImg/cut/%s_2_lf.tif", id));
    imwrite(HForiginal2 * 10, sprintf("outputImg/cut/%s_2_hf.tif", id));

    % part3, LF image activity measure
    % part3.1, raw SML result
    SMLo1 = SML(LForiginal1, SML_N, SML_T, SML_s);
    SMLo2 = SML(LForiginal2, SML_N, SML_T, SML_s);
    % part3.2, softmaxed probability matrices
    [LFPo1, LFPo2] = SoftMax(SMLo1, SMLo2);
    imwrite(LFPo1, sprintf("outputImg/lf/%s_1_pl.tif", id));
    imwrite(LFPo2, sprintf("outputImg/lf/%s_2_pl.tif", id));

    % part4, HF image activity measure
    HFDo1_tmp = vecnorm(HFF(HForiginal1, 52), 2, 3);
    HFDo2_tmp = vecnorm(HFF(HForiginal2, 52), 2, 3);
    [HFDo1, HFDo2] = SoftMax(HFDo1_tmp, HFDo2_tmp);
    imwrite(HFDo1, sprintf("outputImg/hf/%s_1_ph.tif", id));
    imwrite(HFDo2, sprintf("outputImg/hf/%s_2_ph.tif", id));

    % part5, ED
    Ed1 = -log(LFPo1 .* HFDo1);
    Ed2 = -log(LFPo2 .* HFDo2);
    imwrite(Ed1, sprintf("outputImg/ed/%s_1_ed.tif", id));
    imwrite(Ed2, sprintf("outputImg/ed/%s_2_ed.tif", id));

    imwrite(Ed1 > Ed2, sprintf("outputImg/map/%s_raw.tif", id));
    map = FusionMap(Ed1, Ed2, original1, original2);
    imwrite(map, sprintf("outputImg/map/%s_fine.tif", id));

    outImg = (original1 .* ~map) + (original2 .* map);
    imwrite(outImg, sprintf("outputImg/final/%s.tif", id));
end
