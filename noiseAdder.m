function noisyImg = noiseAdder(level, id)
    type = "color";
    if(id(1) ~= 'c')
        type = "grayscale";
    end
    image = im2double(imread(sprintf("testImg/%s/%s.tif", type, id)));
    noisyImg = imnoise(image, 'poisson');
    noisyImg = imnoise(noisyImg, 'speckle', 0.01 * level);
    noisyImg = imnoise(noisyImg, 'gaussian', 0, 0.002 * level);
    noisyImg = imnoise(noisyImg, 'salt & pepper', 0.01 * level);
    imwrite(noisyImg, sprintf("testImg/noisy/n_%d%s.tif", level, id(2:end)));
end