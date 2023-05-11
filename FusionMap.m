% This function create fusion map, don't fuse. Fuse should be done outside

function weightMap = FusionMap(Ed1, Ed2, image1, image2)
    assert(isequal(size(Ed1), size(Ed2)), "Size mismatch");
    assert(isequal(size(image1), size(image2)), "Size mismatch");
    assert(length(size(image1)) <= 3, "Size unexpected");
    assert(length(size(Ed1)) == 2, "Size unexpected");
    
    % some preprocessings
    normImg1 = vecnorm(image1, 2, 3);
    normImg2 = vecnorm(image2, 2, 3);
    [X, Y] = ndgrid(1:size(Ed1, 1), 1:size(Ed1, 2));
    
    % first, based on Eds alone. who large who come.
    weightMap = Ed1 > Ed2;
    
    newMap = zeros(size(weightMap));
    iter = 0;
    maxIter = 100;
    lastUpdate = size(Ed1, 1) * size(Ed1, 2);
    Nc = 10; % NOTE this var determines Es calculation distance
    while (~isequal(weightMap, newMap) && iter < maxIter)
        fprintf("Fusion map entering iteration: %d\n", iter);
        % second, smoothing with Es
        if (iter ~= 0)
            weightMap = newMap;
        end
        iter = iter + 1;
        newMap = arrayfun(@(x, y) PixelDecision( ...
            [x, y], Ed1, Ed2, weightMap, normImg1, normImg2, Nc ...
        ), X, Y);
        update = sum((newMap - weightMap) ~= 0, 'all');
        if (update >= lastUpdate)
            disp("model reaching limit, quit");
            break;
        end
        lastUpdate = update;
        fprintf("Iteration update: %d\n", update);
    end
end
