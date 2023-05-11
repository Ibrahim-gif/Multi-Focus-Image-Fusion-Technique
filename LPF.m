% Ideal low pass filter

% param cut is between 0 and 1. 1 means all pass
function outImg = LPF(inImg, cut)
    imgSize = size(inImg);
    
    % filter matrix
    rowLen = imgSize(1) * cut / 2;
    left = ceil((imgSize(1) + 1)/2 - rowLen);
    right = floor((imgSize(1) + 1)/2 + rowLen);
    colLen = imgSize(2) * cut / 2;
    top = ceil((imgSize(2) + 1)/2 - colLen);
    bottom = floor((imgSize(2) + 1)/2 + colLen);
    lpf = zeros(imgSize);
    lpf(left:right, top:bottom, 1:imgSize(3)) = 1;
    % imshow(lpf);
    
    % fft transform
    fftImg = fftshift(fft2(inImg));
    % imshow(fftImg);
    outImg = ifft2(ifftshift(fftImg .* lpf));
    % imshow(outImg);
end
