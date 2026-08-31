%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


% bm = 
% 
%   struct with fields:
% 
%         diam: 0.0120
%           dx: 8.3226e-06
%           fr: 64.0200
%     PropType: 'OUTSIDE_to_INSIDE_'
%           pz: 0.7683
%        Rbeam: 0
%     RbeamInf: 1
%      RefSurf: 'PLANAR'
%      TypeOld: 'INSIDE_'
%           w0: 5.2985e-05
%        w0_pz: 0.7682
%           wf: [1024×1024 double]
%           wl: 1.3000e-06
%         zRay: 0.0068

% function [wfai, samp] = spidersProper2(wlm, np)
function out = spidersProper5movie(wlm, np)
%        [wfai, samp] = simple_prescription(wlm, np)
%
% Outputs:
% wfai = 2D wave front array intensity
% samp = 2D wave front array sampling distance (m)
%
% Required inputs:
% wlm  = wavelength (m)
% np   = number of pixels

% 2005 Feb     jek  created idl routine
% 2017 Apr 21  gmg  Matlab translation
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
saveDir = '\\Mac\Home\Documents\work\NewEARTH\SubaruPathfinder\simu\';
movieDir = '\\Mac\Home\Documents\work\Gemini\GPI\GPI2.0\CAL2\simu\PassataOutputData\PASSATA mag8 1arcsec seeing\'; 
% aoResMapFile = '\\Mac\Home\Documents\work\Gemini\GPI\GPI2.0\CAL2\simu\gpi2aoResMag8s054.fits';
% aoResMapFile = '\\Mac\Home\Documents\work\Gemini\GPI\GPI2.0\CAL2\simu\gpi2AoRes1arcsecSeeingMag8Snapshot.fits';
apoDir = '\\Mac\Home\Documents\work\NewEARTH\SubaruPathfinder\opticalDesign\Apodizer\Data for fabrication\';
fpmDir = '\\Mac\Home\Documents\work\NewEARTH\SubaruPathfinder\opticalDesign\SPIDERS FPM data\';
% Default parameters
bdf  =    0.100d0          ;  % beam diameter fraction
diam =    7.92d0         ;  % telescope diameter (m)
fr   =   13.901d0          ;  % focal ratio
flm  = diam * fr           ;  % focal length (m)

coro = 1; % 1 for coronagraph

cent = 1; % pupils centering flag
fpmFormat = 'map'; %'anal';
coroBand = 'H'; % select the FPM made by KAUST, if map is selected
fpmType = 'TG';
dis = 0; % 1 for display
pow =0.4; % power scale display for images
detector = 0; % flag for computing detector's image 
cdi = 0; % flag for computing MTF 

% cmap = viridis(256);
% cmap = fake_parula(256);
% cmap = magma(256);
% cmap = plasma(256);
%  cmap = fire(256);
% cmap = inferno(256);
%cmap = cmocean('thermal');
% cmap = cmocean('-matter');
% cmap = cmocean('haline');
% cmap = cmocean('solar');
% cmap = cubehelix(4096,2.65,0.68,1.64,0.4); % custom plasma gamma 0.4 
% cmapGam = cubehelix(4096,0.84,0.2,2.66,0.4); % custom hot gamma 0.4 
cmapLin = cubehelix(512,0.84,0.2,2.66,1);
 set(0,'DefaultFigureColormap',cmapLin);
 % set(0,'DefaultFigureColormap',feval('hot'));

 tel = []; % for tel pupil variable
 apod = []; % apodizer mask variable
 fpm = []; %FPM mask variable
movie = importdata([movieDir 'gpiAOres1arcsecSeeing2secMovie.mat']);
nTime = 1000;
%i=floor(1:99.9:nTime);
sccFrames = zeros(np,np,nTime);
%for i=floor(1:99.9:nTime)
for i=1:nTime % one second only
 fprintf('******* ITERATION %d over %d *******\n',i,nTime); 
bm   = prop_begin(diam, wlm, np, bdf);

% Subaru telescope pupil
pupilGuy = 'Mamadou';
switch pupilGuy
    case 'Mamadou'
        % from Mamadou
        bm   = prop_circular_aperture(bm, diam); % Normalize flux to unobcured pupil, Easier to relate to star mag.
        bm   = prop_define_entrance(bm);
        obsc = 0.2900; % needed later
        if isempty(tel) %  we read the FITS the first time only
            [ bm, tel ] = prop_errormap( bm, 'pupilsbr_nPup1200_kpdiam100_kodiam100_kthick100.fits',...
                'AMPLITUDE','SAMPLING', diam/1200);
        else
            bm = prop_multiply(bm,tel); % apply tel pupil already in memory
        end
    case 'Olivier'
        [bm, obsc]= prop_subaruPupilSpiders(bm, diam / 2.0d0); % Subaru pupil with central obscuration and spider
end

% AO resisuals
[bm, mapOpd] = prop_errormap2(bm,squeeze(movie(:,:,i)),0.045/7.73*diam,'unit','nm');

% bm = prop_errormap(bm, aoResMapFile, 'wavefront', 'sampling', 0.045/7.73*diam, 'nm');

 % Scintillation
% [bm, mapAmp] = prop_psd_errormap(bm,0.01,22/diam,11/3,'AMPLITUDE',1,'TPF'); % 1%RMS scintillation


if dis
displane(bm,'Intensity','Subaru Pupil Intensity',1,6)
displane(bm,'Phase','Subaru Pupil Phase',1,6)
end
bm   = prop_propagate(bm, flm,'snm', 'dummy lens'); % 
bm   = prop_lens(bm, flm, 'F/13.9 convergent beam'); % mimic input telecentric beam
bm   = prop_propagate(bm, flm,'snm', 'Bay #4 focus'); % to Bay#4 focus
if dis
displane(bm,'Intensity','Bay #4 focus',pow,6)
end
bm   = prop_propagate(bm, 772.35e-3,'snm', 'OAE1'); % to OAE1
% displane(bm,'Intensity','OAE1',0.4,4)
[bm, loss]   = prop_lens_aperture(bm, 459.7e-3, 3*25.4e-3,'OAE1', 'disp',dis,'zoom',7, 'power', pow); % OAE1
bm   = prop_propagate(bm, 471.4e-3,'snm', 'DM'); % to DM
if dis
displane(bm,'Intensity','DM',1,6)
end
% bm   = prop_zernikes( bm, 5, 5e-9 ); % apply NCPA on DM

dif = prop_get_distancetofocus(bm) %to bare focus
bm   = prop_propagate(bm, dif+367.2e-3,'snm', 'OAE2'); % to OAE2
% displane(bm,'Intensity','OAE2',pow,4)
[bm , loss]  = prop_lens_aperture(bm, 277.8e-3, 38.1e-3, 'OAE2', 'disp',dis,'zoom',5, 'power',pow); % OAE2
bm   = prop_propagate(bm, 377.2e-3,'snm', 'Apodizer'); % to Apodizer

%% APODIZER


dx = prop_get_sampling(bm);
bmRad = prop_get_beamradius(bm);
switch pupilGuy
    case 'Olivier'
        %  Option 1: Squeezed GPI
        % apod = apodizerGPIsqueeze(bm, obsc);
        
        % Option 2: optimized 1D profile:
        pp2 = importdata('spidersApodizerPolynomial.mat');
        apod = polyRevolutionMask(pp2,np,bmRad,obsc,dx);
        if coro
            bm = prop_multiply(bm,apod); % apply Apodizer
        end
    case 'Mamadou'
        % Option 3: Mamadou apodizer
        if isempty(apod)
        apoFITS = 'pupilsbr_nPup1200_pdiam792_odiam230_thick000_Apod_rMask265.fits'; % grey scale
        % apoFITS = [apoDir 'apodizer9_9umPixel_pupilsbr_nPup1200_pdiam784_odiam238_thick000_Apod_rMask264_EDA_v2.fits']; % half-toned
        [ bm, apod ] = prop_errormap( bm, apoFITS,...
            'AMPLITUDE','SAMPLING', 2*bmRad/1200);
        else
             bm = prop_multiply(bm,apod); % apply Apodizer
        end
end


lambda = prop_get_wavelength(bm); % wavelength (m)
% bm = prop_errormap(bm, aoResMapFile, 'wavefront', 'sampling', 0.045/7.88*12e-3, 'nm');
% [bm, map] = prop_psd_errormap(bm,25e-9,22/diam,11/3,'RMS','TPF'); % Residual WFE
% [bm, mapAmp] = prop_psd_errormap(bm,0.01,22/diam,11/3,'AMPLITUDE',1,'TPF'); % 1%RMS scintillation


%  figure
% imagesc(x*dx*1e3,x*dx*1e3,map); axis equal tight; zoom(2)
% title('Residual AO')
% figure
% imagesc(x*dx*1e3,x*dx*1e3,mapAmp); axis equal tight; zoom(2)
% title('1% RMS Amplitude') 
if dis
displane(bm,'Intensity','Apodized Pupil Intensity',1,6)
end
% displane(bm,'Phase','Apodized Pupil Phase',1,6)
dif2 = prop_get_distancetofocus(bm)
bm   = prop_propagate(bm, dif2,'snm', 'FPM'); % to FPM

if dis
displane(bm,'intensity','FPM PSF',pow,6)
end

%% FPM
% out = fpmMaskProper(bm,radius,tilt,ampGaus, parabx, paraby, charge, piston,apex, nLevel, varargin)
% fpm = fpmMaskProper(bm,453.3e-6/2,1.6*atan(1/64)*4.385/4.518,10,0,0 ,0*4,0,0,0, 'disp',1); % in radian ; GPI FPM J band TG 453.3um
% fpm = fpmMaskProper(bm,453.3e-6/2,0*1.6*atan(1/64),0*10,-0.478,-0.478 ,0*4,0*2*pi,12.4192,0, 'disp',1, 'centering',0); % in radian ; GPI FPM J band TG 453.3um

dx = prop_get_sampling(bm);
dxr = prop_get_sampling_radians(bm);
resel = lambda/bm.diam/dxr*dx % resel in meter


switch fpmFormat
    case 'anal' % amalytic mask
        switch fpmType
            case 'TG'
        radJ= 459.62e-6/2;
        radH= 584.63e-6/2;
        fpm = fpmMaskProperFixQuanta(bm,radH,1.6*atan(1/64),7.8,0,0 ,0*4,1.9999*pi,0,16, 'disp',1, 'centering',cent); % in radian ; GPI FPM J band TG 453.3um
            case 'TGV'
               fpm = fpmMaskProperFixQuanta(bm,resel,1.6*atan(1/64),10,0,0 ,4,0*2*pi,0,16, 'disp',1, 'centering',cent); % in radian ; GPI FPM J band TG 453.3um 
        end
        if coro
            bm = prop_add_phase(bm,fpm/2/pi*lambda); % apply FPM
        end
    case 'map'   % FPM from FITS file
        if isempty(fpm)
        switch coroBand
            case {'i','I'}
                fitsName = [fpmDir 'spidersFpmKaust-Sag-Meter-250nmRes-I.fits'];
            case {'j','J'}
                fitsName = [fpmDir 'spidersFpmKaust-Sag-Meter-250nmRes-J.fits'];
            case {'z','Z'}
                fitsName = [ fpmDir 'spidersFpmKaust-Sag-Meter-250nmRes-Z.fits'];
            case {'h','H'}
                fitsName = [ fpmDir 'spidersFpmKaust-Sag-Meter-250nmRes-H.fits'];
        end
        [bm, fpm] = prop_errormap(bm, fitsName, 'mirr', ...
             'sampling', 250e-9, 'rot', 180);
        else
             bm = prop_add_phase(bm,rot90(fpm,2)*2); % apply FPM in memory
        end
         if cent
             dx = prop_get_sampling(bm);
             x = linspace(-np/2,np/2-1,np)*dx;
             [xx, yy] = meshgrid(x,x);
             tilt2D = (-2*pi*1.6*atan(1/64)/lambda)*yy;
             bm = prop_add_phase(bm,-tilt2D/2/2/pi*lambda); % apply tilt
         end
         
        if dis
            figure
            dx = prop_get_sampling(bm);
            x = linspace(-np/2,np/2-1,np)*dx;
            imagesc(x*1e3,x*1e3,fpm); axis square tight; zoom(4)
            colormap(cubehelix(256,0.84,0.2,2.66,1));
            cbh = colorbar;
            cbh.Label.String = 'Phase [rad]';
            
            title('FPM')
            xlabel('[mm]')
            ylabel('[mm]')
        end
end
%%
bm   = prop_propagate(bm, 358.69e-3+5.166e-3, 'snm','Lens 0'); % to lens0 object principal point
% displane(bm,'intensity','Lens 0',0.4,3)
[bm, loss]   = prop_lens_aperture(bm, 285.99e-3, 40e-3,'Lens 0', 'disp',dis,'power',pow, 'zoom',1); % Lens0
bm   = prop_propagate(bm, 382.67e-3 -0.022e-3,'snm','Lyot'); % from lens0 image principal plane
if dis
   displane(bm,'intensity','Lyot pupil before stop',pow,3)%
end
%% Lyot stop
pupDia = 2*2.034e-3; % 4.049e-3; % 7.115e-3 % in meter at Lyot plane
 sep = 1.6*pupDia;
% sep = 3.615e-3;

   % lyot = sccLyotMaskSubaruProper(bm,[1 1/43*24/14]*pupDia, sep, -pi/2, 'obsDia', obsc*pupDia, 'disp',dis, 'lyotMargins', [0.98 1 1.2], 'center', [0, 1.6*pupDia/2], 'power',pow, 'zoom',3);
   lyot = sccLyotMaskSubaruProperGaussPin(bm,[1 1/20]*pupDia, sep, -pi/2, 'obsDia', obsc*pupDia, 'disp',dis, ...
       'lyotMargins', [0.95 1.15 2], 'center', [0, cent*1.6*pupDia/2], 'power',pow, 'zoom',3);

  fprintf('Pinhole Transmission    = %8.5e\n',lyot.pinholeTrans);
  fprintf('Main Pupil Transmission = %8.5e\n',lyot.mainPupilTrans);


if coro
bm = prop_multiply(bm, lyot.mask); % apply Lyot stop
end
% phLens = standardZemaxSurface(np,dx,pupDia/2,0.05,0,[0. -1.6*pupDia]); % OPD in meter
% bm = prop_add_phase(bm,phLens);

% figure
% imagesc(x*dx*1e3,x*dx*1e3,100*phLens+abs(pupima).^2); axis equal tight; zoom(2)
% title('Lyot pupil+Lens')
%%

bm   = prop_propagate(bm, 75e-3,'snm','Chopper'); %
if dis
displane(bm,'intensity','Chopper plane',pow,3)
end
bm   = prop_propagate(bm, 113.16e-3-0.022e-3,'snm','Lens 1'); %
[bm, loss]   = prop_lens_aperture(bm, 285.99e-3, 40e-3,'Lens 1', 'disp',dis, 'zoom',2, 'power', pow); % Lens1
dif = prop_get_distancetofocus(bm) %to bare focus
bm   = prop_propagate(bm, dif,'snm','SCC Focus'); %
if dis
displane(bm,'intensity','SCC Focus',pow,dis)
% displane(bm,'intensity','SCC Focus',1,1)
end

% compute contrast profile
if dis
cont = prop_psfContrast(bm,coro);
end

if detector
    % simulate C-RED2 detector
    pixelSize = 15e-6;
    nPix = 512;
    gain = 33e3/2^14; % e-/ADU
    QE = 0.7;
    ron = 25;
    expo = 1/1; % seconds
    mag = 6;
    
    camIma =  prop_detector(bm, pixelSize, nPix, gain,QE, ron, expo ,mag, obsc, 'disp',1, 'snm', 'C-RED2');
end

if coro
    if cdi
       mtf = cdiOTF(bm,1*.9,1.575,pi/2,'disp',dis);
    end
else
    cont.psf2D = (prop_get_amplitude(bm)).^2; % save non-coronagraphic PSF and profile for reference
    save('nonCoronPSF.m','cont')
end

%% final sampling
finalSampling =  1.1513e-05;

dx = prop_get_sampling(bm)
 mag = dx/ finalSampling;
 field = prop_magnify(prop_get_wavefront(bm), mag, 'size_out', 1024, 'conserve');
 sccFrames(:,:,i) = abs(field).^2;

% Calculate 2D wave front array intensity and sampling distance (m)
 % [wfai, samp] = prop_end(bm);

end
out = sccFrames;
waveStr = num2str(wlm*1e9,'%d');
% save SSC Frame cube
save([saveDir 'spidersSccFrames_Mag8_Seeing1_Fps1kHz_H_band_FPM_' waveStr 'nm.mat'],'sccFrames','-v7.3');

end                     % function simple_prescription
