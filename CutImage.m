function [LFImg1, HFImg1, LFImg2, HFImg2] = CutImage(inputImg1, inputImg2)
    percentageLow = 0.01;
    percentageHigh = 0.015;
    top = 1;
    bottom = 0;
    totalPixel = sum(inputImg1, 'all') + sum(inputImg2, 'all');
    cut = 0.5;
    while(1)
        cut = (top + bottom) / 2;
        HFImg1 = abs(HPF(inputImg1, cut));
        HFImg2 = abs(HPF(inputImg2, cut));
        HFPixel = sum(HFImg1, 'all') + sum(HFImg2, 'all');
        HFPercentage = HFPixel / totalPixel;
        if (HFPercentage < percentageLow)
            % cut should be lower
            top = cut;
        elseif (HFPercentage > percentageHigh)
            % cut should be higher
            bottom = cut;
        else
            break;
        end
    end
    LFImg1 = abs(LPF(inputImg1, cut));
    LFImg2 = abs(LPF(inputImg2, cut));
end