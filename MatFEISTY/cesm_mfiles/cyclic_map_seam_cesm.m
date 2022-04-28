function [lat2,lon2,gvar2] = cyclic_map_seam_cesm(lat,lon,gvar)

lon2 = lon;
lat2 = lat;
gvar2 = gvar;

% Determine which direction is latitude
[ni,nj] = size(lon);
% if (ni > nj)
  lon2(ni+1,:) = lon(ni,:)+1;
  lat2(ni+1,:) = lat(ni,:);
  gvar2(ni+1,:) = gvar(ni,:)+eps;
% else
%   lon2(:,nj+1) = lon(:,nj)+1;
%   lat2(:,nj+1) = lat(:,nj);
%   gvar2(:,nj+1) = gvar(:,nj)+eps;
% end
end
