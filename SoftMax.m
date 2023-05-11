% Compute the soft max between two matrices

function [out1, out2] = SoftMax(in1, in2)
    out1 = exp(in1) ./ (exp(in1) + exp(in2));
    out2 = exp(in2) ./ (exp(in1) + exp(in2));
    % damn its simple
end
