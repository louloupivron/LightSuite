function cordvol = cordDownsampleVolume(cordvol, opts)
%CORDDOWNSAMPLEVOLUME Resample cord volume to the registration grid if needed.
%   Skips imresize3 when cordvol is already at the target size (e.g. after
%   streaming load in readSpinalCordSample).
if isfield(opts, 'orisize') && ~isempty(opts.orisize)
    nativeSize = opts.orisize(:).';
else
    nativeSize = size(cordvol, 1:3);
end
[resfac, targetSize] = cordVolumeDownsampleSpec(nativeSize, ...
    opts.sampleres, opts.registrationres);
if isequal(size(cordvol, 1:3), targetSize)
    return
end
Nchan = size(cordvol, 4);
finvoluse = zeros([targetSize Nchan], 'uint16', 'like', cordvol);
for ii = 1:Nchan
    finvoluse(:, :, :, ii) = imresize3(cordvol(:, :, :, ii), 'Scale', resfac);
end
cordvol = finvoluse;
end
