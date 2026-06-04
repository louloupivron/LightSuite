function [ny, nx, nz, planesInTime, preferImread] = tiffChannelStackDims(path)
%TIFFCHANNELSTACKDIMS Height, width, and depth for one channel TIFF stack.
%   preferImread is true when numel(imfinfo)>1 (classic multi-IFD stack).
%   When imfinfo reports a single directory (OME-TIFF, SubIFDs, etc.),
%   use Bio-Formats sizeZ/sizeT instead and set preferImread false.
%   planesInTime is true when stack depth is stored along Bio-Formats T.
tinfo  = imfinfo(path);
nPages = numel(tinfo);
ny     = tinfo(1).Height;
nx     = tinfo(1).Width;
if nPages > 1
    nz             = nPages;
    planesInTime   = false;
    preferImread   = true;
    return;
end
preferImread = false;
datainfo = BioformatsImage(path);
ny       = datainfo.height;
nx       = datainfo.width;
sizeZ    = datainfo.sizeZ;
if isprop(datainfo, 'sizeT')
    sizeT = datainfo.sizeT;
else
    sizeT = 1;
end
if sizeZ == 1 && sizeT > 1
    nz           = sizeT;
    planesInTime = true;
else
    nz           = sizeZ;
    planesInTime = false;
end
end
