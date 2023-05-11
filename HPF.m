% Ideal high pass filter

% param cut is between 0 and 1. 1 means no pass
function outImg = HPF(inImg, cut)
    imgSize = size(inImg);
    
    % filter matrix
    rowLen = imgSize(1) * cut / 2;
    left = ceil((imgSize(1) + 1)/2 - rowLen);
    right = floor((imgSize(1) + 1)/2 + rowLen);
    colLen = imgSize(2) * cut / 2;
    top = ceil((imgSize(2) + 1)/2 - colLen);
    bottom = floor((imgSize(2) + 1)/2 + colLen);
    hpf = ones(imgSize);
    hpf(left:right, top:bottom, 1:imgSize(3)) = 0;
    % imshow(lpf);
    
    % fft transform
    fftImg = fftshift(fft2(inImg));
    % imshow(fftImg);
    outImg = ifft2(ifftshift(fftImg .* hpf));
    % imshow(outImg);
end