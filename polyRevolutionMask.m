%% Compute a centro-symmetric mask defined by a polynimial radial profile
% np : np x np pixels map
% p : vector with polynomial coefficient
% dx sampling in meter
% pupRad telescope pupil semi-diameter (m)
% obsc : central obscuration ratio

function out = polyRevolutionMask(pol, np,pupRad,obsc,dx)
%maxDegree = 7;
x = linspace(-np/2,np/2-1,np)*dx; % in meter
[xx, yy] = meshgrid(x,x);
r = sqrt(xx.^2+yy.^2)/pupRad; % normalized radius

%p = zeros(1,maxDegree+1);
%p(1:numel(pol)) = pol;  % just in case pol is too short

prof = 0;
for i=1:numel(pol)
  prof = prof + pol(i)*r.^(i-1);
  % prof = p(1)*r.^7 + p(2)*r.^6 + p(3)*r.^5 + p(4)*r.^4 + p(5)*r.^3 + p(6)*r.^2 + p(7)*r + p(8) ;
end  
 % mask = prof/max(prof(:)) .* prop_ellipse(bm,pupRad,pupRad) .* prop_ellipse(bm,obsRad*pupRad,obsRad*pupRad, 'dark');

   % prof(r>1)=0;
   % prof(r<obsc) = 0;
   mask = prof/max(prof(r>=obsc & r<=1));
   mask(mask<0) = 0;
out = mask;
end