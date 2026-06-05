function [resfac, targetSize] = cordVolumeDownsampleSpec(nativeSize, sampleres, registrationres)
%CORDVOLUMEDOWNSAMPLESPEC Pixel scale and target size for cord registration grid.
%   [resfac, targetSize] = cordVolumeDownsampleSpec(nativeSize, sampleres, registrationres)
%   nativeSize is [dim1 dim2 dim3] in pixels; sampleres and registrationres are
%   voxel sizes in micrometers (scalar or [x y z]).
nativeSize = nativeSize(:).';
sampleres = localNormalizeRes(sampleres);
registrationres = localNormalizeRes(registrationres);
resfac = sampleres ./ registrationres;
targetSize = max(1, ceil(nativeSize .* resfac));
end

%--------------------------------------------------------------------------
function r = localNormalizeRes(r)
r = r(:).';
if numel(r) == 1
    r = repmat(r, 1, 3);
end
end
