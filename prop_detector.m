function out = prop_detector(bm, pixelSize, nPix, gain,QE, ron, expo ,mag, obsc, varargin)
% Scale image in electrons based on star mag, expo and QE(assume bm.diam for telescope
% diameter)
% Resample image, rescale to keep sum cte. 
% Crop to nPix (either nx or ny, your choice) pixels
% Add photon noise and read-out noise
% Quantize image on 2^nBit level

p = inputParser;
p.addParameter('disp',0, @isnumeric);
p.addParameter('snm','Detector', @ischar);
p.parse(varargin{:});

% f0 is Zero Point in ph/s/m2
f0 = [13.5 5.63 4.38 3.22 2.82 1.51] * 1e9; % Zero Point for R,Z',Y,J,H,K
area = pi*(bm.diam/2)^2; % entrance is defined before obscuration
% nb of photon-events/s at entrance pupil (sum=1)
nph = f0*area*10.^(-0.4*mag)*expo*QE;

lambda = prop_get_wavelength(bm);
dx = prop_get_sampling(bm);
dxr = prop_get_sampling_radians(bm);
psf = (prop_get_amplitude(bm)).^2;

if lambda<0.73e-6 % R
    nel0 = nph(1);
    band = 'R';
elseif lambda<0.98e-6 % I or Z'
    nel0 = nph(2);
    band = 'Z';
elseif lambda<1.1e-6 % Y
    nel0 = nph(3);
    band = 'Y';
elseif lambda<1.42e-6 % J
    nel0 = nph(4);
    band = 'J';
elseif lambda<1.9e-6 % H
    nel0 = nph(5);
    band = 'H';
elseif lambda<2.5e-6
    nel0 = nph(6);
    band = 'K';
else
    disp('Incorrect wavelenght! Assuming J band.')
    nel0 = nph(4);
    band = 'J';
end
      
% remaining nb of photo-events after corono and after optics throughput

tauOpt = 0.45*10/17.45; % T=0.45 for Tel*AO188*BSwitcher*SPIDERS optics, 10% band pass in J
% scale psf in photo-events:
psf2 = psf*nel0*tauOpt;

% sample image on pixels
psfSamp = prop_magnify( psf2, dx/pixelSize,'CONSERVE', 'SIZE_OUT', nPix );

% add photon and read-out noise
psfNoise = 1e12*imnoise(psfSamp*1e-12,'poisson') + ron*randn(size(psfSamp));

% scale image to ADU and quantize it (no saturation)
out = floor(psfNoise/gain); % gain is in e-/ADU

if p.Results.disp
figure
x = linspace(-nPix/2,nPix/2-1,nPix);
samp = dxr/lambda*bm.diam/dx*pixelSize;
imagesc(x*samp,1:nPix,out); axis square tight;
colormap(cubehelix(256,0.84,0.2,2.66,1));
hcb = colorbar;
hcb.Label.String = 'Intensity [ADU]';
xlabel('[\lambda/D]')
ylabel('[detector pixels]')
title([sprintf('%s Image | %3.2f pixels/resel @ %3.2f',p.Results.snm, 1/samp, lambda*1e6), '\mum', ...
    sprintf(' | %s = %2.1f | Expo. = %3.2f s',band,mag,expo)])
end
end



