% This one is the sub method for FusionMap,
% for energy comparison based on pixel's decision.

function choice = PixelDecision(p, Ed1, Ed2, cMap, img1, img2, N)
    p1Norm = img1(p(1), p(2));
    p2Norm = img2(p(1), p(2));
    xWindow = max(1, p(1)-N):min(p(1)+N, size(Ed1, 1));
    xCenter = min(p(1), 1 + N);
    yWindow = max(1, p(2)-N):min(p(2)+N, size(Ed1, 2));
    yCenter = min(p(2), 1 + N);
    windowChoose = cMap(xWindow, yWindow);
    windowImg1Norm = img1(xWindow, yWindow);
    windowImg2Norm = img2(xWindow, yWindow);
    windowFg = sqrt((p1Norm - windowImg2Norm).^2 + ...
                    (p2Norm - windowImg1Norm).^2);
    windowC = sqrt((p1Norm - windowImg1Norm).^2 + ...
                   (p2Norm - windowImg2Norm).^2);
    windowDMap = zeros(size(windowC));
    windowDMap(xCenter, yCenter) = 1;
    windowDMap = bwdist(windowDMap);
    c = 0.0001;
    Es = (1 ./ windowDMap) .* (windowFg + c) ./ (windowC + c);
    Es(xCenter, yCenter) = 0;
    totalE1 = Ed1(p(1), p(2)) + sum(Es .* windowChoose, 'all');
    totalE2 = Ed2(p(1), p(2)) + sum(Es .* ~windowChoose, 'all');
    choice = totalE1 > totalE2;
end
