% SML: Sum of Modified Laplacian
% Low freq image activity measure. This function does not go to softmax.

% check paper for what they stand, if not obvious
% paper didn't say that, but I will use zero padding.
function smlMat = SML(inImg, N, T, s)
    imgSize = size(inImg);
    padImg = padarray(inImg, [s, s], 0, 'both');
    % step 1: ML matrix.
    mlMat = zeros(imgSize(1), imgSize(2));
    % per-pixel base calculation
    for x = 1:imgSize(1)
        for y = 1:imgSize(2)
            rowWise = norm(reshape(2*padImg(x+s, y+s, :) ...
                                   - padImg(x+s-s, y+s, :) ...
                                   - padImg(x+s+s, y+s, :), 1, []));
            colWise = norm(reshape(2*padImg(x+s, y+s, :) ...
                                   - padImg(x+s, y+s-s, :) ...
                                   - padImg(x+s, y+s+s, :), 1, []));
            mlMat(x, y) = rowWise + colWise;
        end
    end
    % imshow(mlMat);
    padMl = padarray(mlMat, [N, N], 0, 'both');
    smlMat = zeros(imgSize(1), imgSize(2));
    % again, per-pixel, Pad mlMat with zeros again
    for x = 1:imgSize(1)
        for y = 1:imgSize(2)
            mls = padMl(x:x+2*N, y:y+2*N);
            mls(mls < T) = 0;
            smlMat(x, y) = sum(mls, 'all');
        end
    end
    % imshow(smlMat);
end
