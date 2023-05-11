% NOTE: add the "objective_metrics" folder to your path.
function m = metrics(id)
    type = "color";
    if (id(1) == 'g')
        type = "grayscale";
    elseif (id(1) == 'n')
        type = "noisy";
    end
    sourceId = id;
    if (id(1) == 'n')
        sourceId(1) = 'c';
    end
    in1 = im2gray(imread(sprintf('testImg/%s/%s_1.tif', type, sourceId)));
    in2 = im2gray(imread(sprintf('testImg/%s/%s_2.tif', type, sourceId)));
    fuse = im2gray(imread(sprintf('outputImg/final/%s.tif', id)));
    
    % Ref: Comments on "Information measure for performance of image fusion"
    % By M. Hossny et al.
    % Electronics Letters Vol. 44, No.18, 2008
    MI = metricMI(in1, in2, fuse, 1);
    
    % Slightly modified by Yu Liu after a careful check according to original publication [Oct 2019]  
    % Ref: Objective Image Fusion Performance Measure, Electronics Letters (2000)
    % Ref: Objective Pixel-level Image Fusion Performance Measure, Proc. SPIE 4051, 89 (2000)
    % by C. Xydeas 
    QG = metricXydeas(in1, in2, fuse);
    
    % Slightly modified by Yu Liu after a careful check according to original publication [Oct 2019]

    % Ref: A novel similarity based quality metric for image fusion, Information Fusion, Vol.9, pp156-160, 2008
    % by Cui Yang et al.
    QY = metricYang(in1, in2, fuse);
    
    % Ref: A new automated quality assessment algorithm for image fusion, Image and Vision Computing, 27 (2009) 1421-1432 
    % By Yin Chen et al.
    % 
    QC = metricChenBlum(in1, in2, fuse);
    
    m = [MI, QG, QY, QC];
end