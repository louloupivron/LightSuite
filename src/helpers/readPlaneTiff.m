function img = readPlaneTiff(path)
%READPLANETIFF Read a single 2D plane from a one-page TIFF.
%   img = readPlaneTiff(path) tries, in order: Tiff class, imread, bfopen.
%   Returns a 2D numeric array (typically uint16).

path = char(path);
errors = {};

try
    tf = Tiff(path, 'r');
    img = read(tf);
    close(tf);
    img = localTo2dPlane(img);
    return
catch ME
    errors{end+1} = sprintf('Tiff: %s', ME.message); %#ok<AGROW>
end

try
    img = imread(path);
    img = localTo2dPlane(img);
    return
catch ME
    errors{end+1} = sprintf('imread: %s', ME.message); %#ok<AGROW>
end

try
    data = bfopen(path);
    img = data{1}{1, 1};
    if iscell(img)
        img = img{1};
    end
    img = localTo2dPlane(img);
    return
catch ME
    errors{end+1} = sprintf('bfopen: %s', ME.message); %#ok<AGROW>
end

error('readPlaneTiff:ReadFailed', ...
    'Could not read plane TIFF:\n  %s\n%s', path, strjoin(errors, newline));
end

%--------------------------------------------------------------------------
function img = localTo2dPlane(img)
if ndims(img) > 2
    img = img(:, :, 1);
end
if ~isa(img, 'uint16')
    img = im2uint16(img);
end
end
