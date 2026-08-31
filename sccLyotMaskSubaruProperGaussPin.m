function out = sccLyotMaskSubaruProperGaussPin(bm,dias, sep, pa, varargin)
% Creates a double hole Lyot mask for SCC (using TG or TGV focal plane masks)
% Does NOT apply the mask to the wavefront !!!!!!!!!!!

% INPUTS:
%   bm = current PROPER beam
%   dias = 2-element array with holes diameters in meter
%   sep = separation of the 2 holes in meters
%   pa = position angle of the pinhole in radians
%   (opttional) obsDia = Lyot mask central obscuration diameter in meter (none by defualt)
%   (optional) disp = 1 for display, 0 for none (defualt)

% OUTPUTS:
% out.mask = lyot mask image
% out.mainPupiTrans = transmission of the main pupil (leakage)
% out.pinholeTrans = transmission of the pinhole

p = inputParser;
p.addParameter('obsDia',0, @isnumeric);
p.addParameter('disp',0, @isnumeric);
p.addParameter('center',[0 0], @isnumeric);
p.addParameter('lyotMargins',[1 1 1], @isnumeric);
p.addParameter('power',0.4, @isnumeric);
p.addParameter('zoom',3, @isnumeric);
p.parse(varargin{:});
obsDia = p.Results.obsDia;
cc = p.Results.center;
lyotMargins = p.Results.lyotMargins;
pow = p.Results.power;
zoo = p.Results.zoom;

pupRad = dias(1)/2;
% main pupil
mainPup = prop_ellipse(bm,lyotMargins(1)*dias(1)/2,lyotMargins(1)*dias(1)/2, 'xc', cc(1), 'yc', cc(2));
if obsDia % if central obscuration
    mainPup = mainPup .* prop_ellipse(bm,lyotMargins(2)*obsDia/2,lyotMargins(2)*obsDia/2,'xc',cc(1),'yc',cc(2),'dark');
end

% spiders
a = tand(51.5); % angle of the spiders
b2 = -30.7/273*pupRad*2;
b1 = -22.3/273*pupRad*2;
thic = abs(b1-b2); % spider thickness in meter

sp1x = [0 0 1 1]*pupRad;
sp1y = [b2 b1 a*pupRad+b1 a*pupRad+b2] + 0.5*thic*(lyotMargins(3)-1)*[-1 1 1 -1]; 
sp2y = -sp1y;

sp3x = -sp1x;
sp3y = -sp1y;
sp4y = -sp2y;

spider1 = prop_irregular_polygon( bm, sp1x+cc(1), sp1y+cc(2), 'DARK');
spider2 = prop_irregular_polygon( bm, sp1x+cc(1), sp2y+cc(2), 'DARK');
spider3 = prop_irregular_polygon( bm, sp3x+cc(1), sp3y+cc(2), 'DARK');
spider4 = prop_irregular_polygon( bm, sp3x+cc(1), sp4y+cc(2), 'DARK');

mainPup = mainPup.*spider1.*spider2.*spider3.*spider4;

% pinhole
phx = sep*cos(pa)+cc(1);
phy = sep*sin(pa)+cc(2);

% phProfile = 'gauss';
 phProfile = 'top hat';
% phProfile = 'snowman';
% phProfile = 'dart';
switch phProfile
    case 'top hat'
phole = prop_ellipse(bm,dias(2)/2,dias(2)/2,'xc',phx,'yc',phy);

    case 'snowman'
        dias3 = dias(2)*20/14;
phole = prop_ellipse(bm,dias(2)/2,dias(2)/2,'xc',phx+dias(2)/2,'yc',phy)...
    + (14/20)^2*prop_ellipse(bm,dias3/2,dias3/2,'xc',phx-dias3/2,'yc',phy);
 case 'dart'
     ratio = 20/10; %20/14;   % 2nd pinhole is bigger
     dias3 = dias(2)*ratio;
phole = (1-1/ratio^2)* prop_ellipse(bm,dias(2)/2,dias(2)/2,'xc',phx,'yc',phy)...
    + 1/ratio^2 * prop_ellipse(bm,dias3/2,dias3/2,'xc',phx,'yc',phy);
    case 'gauss'
dim = prop_get_gridsize();
dx = prop_get_sampling(bm);
% dxr = prop_get_sampling_radians(bm);
% lambda = prop_get_wavelength(bm); % wavelength (m)
% resel = lambda/bm.diam/dxr*dx; % resel in meter
x = linspace(-dim/2,dim/2-1,dim)*dx; %  matches PROPER origin in meter
[xx, yy] = meshgrid(x-phx,x-phy);
rho = sqrt(xx.^2+yy.^2);
% phole = 1*exp(-0.5*(rho/(dias(2)/2)).^2);
rho2 = rho/(dias(2)/2)*1.61; % normalize rho to FWHM
phole = (2*besselj(1,rho2)./rho2).^2; % Airy disc pinhole
end
out.mask = 1*mainPup + 1*phole;

% compute transmission
pupima = (prop_get_amplitude(bm)).^2;
out.mainPupilTrans = sum(sum(mainPup.*pupima))/sum(pupima(:));
out.pinholeTrans = sum(sum(phole.*pupima))/sum(pupima(:));

if p.Results.disp
    np = size(pupima,1);
    x = linspace(-np/2,np/2-1,np);
    dx = prop_get_sampling(bm);
    cmap = cubehelix(4096,0.84,0.2,2.66,pow);
    figure
    imagesc(x*dx*1e3,x*dx*1e3,((1-out.mask).*pupima).^1); axis equal tight; zoom(zoo)
    colormap(cmap); colorbar;
    title(sprintf('Reflected Lyot pupil | Pupil Trans = %0.3e | Pinhole Trans = %0.3e', out.mainPupilTrans, out.pinholeTrans))
    xlabel('[mm]')
    ylabel('[mm]')
    figure
    imagesc(x*dx*1e3,x*dx*1e3,(out.mask.*pupima).^1); axis equal tight; zoom(zoo)
    colormap(cmap); colorbar;
    title(sprintf('Transmitted Lyot pupil | Pupil Trans = %0.3e | Pinhole Trans = %0.3e', out.mainPupilTrans, out.pinholeTrans))
    xlabel('[mm]')
    ylabel('[mm]')
end
end

