% Pyramid Tilt Parabola Vortex FPM

function out = fpmMaskProperFixQuanta(bm,radius,tilt,ampGaus, parabx, paraby, charge, piston,apex, nLevel, varargin)
% output the phase in radians
% radius of central part in meter
% tilt slope is in rad
% parabolic sag in radian of phase in x and y 
% charge of the vortex
% nLevel = nb of levels for quantization. 0 means no quantization.
% apex pyr angle is actually the PV phase
% 'disp',1 for display (optional) 


p = inputParser;
p.addParameter('disp',0, @isnumeric);
p.addParameter('centering',0, @isnumeric);
p.parse(varargin{:});

%x = linspace(-(dim-1)/2,(dim-1)/2,dim)-.5; % -.5 to match PROPER origin
dim = prop_get_gridsize();
dx = prop_get_sampling(bm);
dxr = prop_get_sampling_radians(bm);
lambda = prop_get_wavelength(bm); % wavelength (m)
resel = lambda/bm.diam/dxr*dx; % resel in meter
x = linspace(-dim/2,dim/2-1,dim)*dx; %  matches PROPER origin in meter
[xx, yy] = meshgrid(x,x);
rho = sqrt(xx.^2+yy.^2);
theta = atan2(-xx,-yy);

coreMask = prop_ellipse(bm,radius,radius);
vortexMask = prop_ellipse(bm,radius,radius, 'DARK');

circCenter = rho<radius;
squaCenter = abs(xx)<radius & abs(yy)<radius;
center = circCenter;

tilt2D = (-2*pi*tilt/lambda)*yy;
%tilt2D = -tilt/radius*(xx+yy)/sqrt(2); % ddiagonal tilt

% Pyramid
pyr = zeros(dim,dim);
nFace = 3;
switch nFace
    case 4
        pyr1 = apex/resel*xx;
        pyr3 = apex/resel*yy;
        pyr((yy>=-xx) & (yy<xx)) = -pyr1((yy>=-xx) & (yy<xx)) ;
        pyr((yy<-xx) & (yy>=xx))=  pyr1((yy<-xx) & (yy>=xx));
        pyr(yy>=abs(xx))= -pyr3(yy>=abs(xx));
        pyr(yy<-abs(xx))= pyr3(yy<-abs(xx)) ;
    case 3
        t1 = 0;
        t2 = pi/4; % 2*pi/3;
        t3 = -pi/4; % -2*pi/3;
        pyr1 = (2*pi*apex/lambda)*rho.*cos(theta-t1);
        pyr2 = (2*pi*apex/lambda)*rho.*cos(theta-t2);
        pyr3 = (2*pi*apex/lambda)*rho.*cos(theta-t3);
        pyr(theta>=-pi/3 & theta<pi/3) = pyr1(theta>=-pi/3 & theta<pi/3);
        pyr(theta>=pi/3 & theta<pi) = pyr2(theta>=pi/3 & theta<pi);
        pyr(theta>=-pi & theta<-pi/3) = pyr3(theta>=-pi & theta<-pi/3);
end
% pyr = pyr - min(pyr(rho<radius)); % to avoid negative sag

% sy = radius*5; % radius = 1 lamb/D
% sx = 0.7*sy;
%  gaus = ampGaus*exp(-0.5*(xx.^2/sx^2+yy.^2/sy^2));
 % s = resel*5; % TGV
 s = resel*2; % TG
gaus = ampGaus*exp(-0.5*(rho/s).^2);
% gaus = gaus - min(gaus(center));

parab = -parabx/radius^2*xx.^2 - paraby/radius^2*yy.^2;
% parab = parab - min(parab(center));

% if charge == 0
%     piston = 0; % no piston if no vortex
% end
vortex = charge*(theta)+piston; 

% vortex(center)= tilt2D(center) + parab(center) + gaus(center) + pyr(center);

core = (tilt2D + parab + gaus + pyr);
core = coreMask.*(core - min(core(coreMask>0))); % avoid negative sag
out = vortexMask.*vortex + core; 
if nLevel>0 % quantization
    refLambda = 1.590e-6;
   hStep = refLambda/lambda*2*pi/nLevel;
   out = hStep*floor(mod(out,2*pi)/hStep); 
end

if p.Results.disp % display before centering
    figure
    imagesc(x*dxr/lambda*bm.diam/dx,x*1e3,out); axis square tight; zoom(4)
    colormap(cubehelix(256,0.84,0.2,2.66,1));
    cbh = colorbar;
    cbh.Label.String = 'Phase [rad]';
    
    title('FPM')
    xlabel('[\lambda/D]')
    ylabel('[mm]')
end
 out = out - p.Results.centering*tilt2D/2; % to center both pupils
 
end