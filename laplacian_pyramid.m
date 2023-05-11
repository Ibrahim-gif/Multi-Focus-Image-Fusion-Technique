close all
%imagesfiles = dir(fullfile('MFIF-master\sourceimages\color\', '*.tif'));
%nfiles = length(imagesfiles);
%for i=1:nfiles
%    currentfilename = fullfile('MFIF-master\sourceimages\color\', imagesfiles(i).name);
%    currentimage = imread(currentfilename);
%    figure
%    image(currentimage);
%end
A = imread("MFIF-master\sourceimages\color\c_01_1.tif");
figure
image(A)
figure
A1 = imread("MFIF-master\sourceimages\color\c_01_2.tif");
image(A1)
A2 = imread("MFIF-master\Results\color\c_01_brw.tif");
figure
image(A2)
M = size(A, 1);
N = size(A, 2);
[M, N];

pyramid = multiresolutionPyramid(A);
visualizePyramid(pyramid);

function mrp = multiresolutionPyramid(A,num_levels)
%multiresolutionPyramid(A,numlevels)
%   mrp = multiresolutionPyramid(A,numlevels) returns a multiresolution
%   pyramd from the input image, A. The output, mrp, is a 1-by-numlevels
%   cell array. The first element of mrp, mrp{1}, is the input image.
%
%   If numlevels is not specified, then it is automatically computed to
%   keep the smallest level in the pyramid at least 32-by-32.

%   Steve Eddins
%   MathWorks

A = im2double(A);

M = size(A,1);
N = size(A,2);

if nargin < 2
    lower_limit = 32;
    num_levels = min(floor(log2([M N]) - log2(lower_limit))) + 1;
else
    num_levels = min(num_levels, min(floor(log2([M N]))) + 2);
end

mrp = cell(1,num_levels);

smallest_size = [M N] / 2^(num_levels - 1);
smallest_size = ceil(smallest_size);
padded_size = smallest_size * 2^(num_levels - 1);

Ap = padarray(A,padded_size - [M N],'replicate','post');

mrp{1} = Ap;
for k = 2:num_levels
    mrp{k} = imresize(mrp{k-1},0.5,'lanczos3');
end

mrp{1} = A;
end
function tiles_out = visualizePyramid(p)
% Steve Eddins
% MathWorks

M = size(p{1},1);
N = size(p{1},2);

for k = 1:numel(p)
    Mk = size(p{k},1);
    Nk = size(p{k},2);
    Mpad1 = ceil((M - Mk)/2);
    Mpad2 = M - Mk - Mpad1;
    Npad1 = ceil((N - Nk)/2);
    Npad2 = N - Nk - Npad1;

    A = p{k};
    A = padarray(A,[Mpad1 Npad1],0.5,'pre');
    A = padarray(A,[Mpad2 Npad2],0.5,'post');
    p{k} = A;
end

tiles = imtile(p,'GridSize',[NaN 2],'BorderSize',20,'BackgroundColor',[0.3 0.3 0.3]);
imshow(tiles)

if nargout > 0
    tiles_out = tiles;
end
end
function f = lanczos3(x)
% See Graphics Gems, Andrew S. Glasser (ed), Morgan Kaufman, 1990,
% pp. 157-158.

f = (sin(pi*x) .* sin(pi*x/3) + eps) ./ ((pi^2 * x.^2 / 3) + eps);
f = f .* (abs(x) < 3);
end