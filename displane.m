function displane(bm,what,titleStr, pow,zoo)
np = prop_get_gridsize();
x = linspace(-np/2,np/2-1,np);

dx = prop_get_sampling(bm);
dxr = prop_get_sampling_radians(bm);
lamb = prop_get_wavelength(bm);
dfoc = abs(prop_get_distancetofocus(bm));

switch what
    case {'phase', 'Phase'}
        wf = prop_get_wavefront(bm);
        intens = (abs(wf).^2);
        ph = prop_get_phase(bm); % in rad
        ima = ph;
        ima(intens<max(intens(:))/10000)=0; % do not diplay phase where there is no light
        cmap = cubehelix(256,0.84,0.2,2.66,1);
    otherwise
        wf = prop_get_wavefront(bm);
        if isnumeric(pow)
           ima = (abs(wf).^2);
           if pow == 1
               cmap = cubehelix(512,0.84,0.2,2.66,1);
           else
               cmap = cubehelix(4096,0.84,0.2,2.66,pow); % more colors because of gamma compression/expansion
           end
        else
            ima = log(abs(wf).^2);
            cmap = cubehelix(512,0.84,0.2,2.66,1);
        end
end
figure
if (bm.w0<=0.5e-3) && (dfoc<10e-3) % Focal plane
  imagesc(x*dxr/lamb*bm.diam,x*dx*1e3,ima); axis square tight; zoom(zoo)
  xlabel('[\lambda/D]')
else % pupil plane (or others)
    imagesc(x*dx*1e3,x*dx*1e3,ima); axis equal tight; zoom(zoo)
    xlabel('[mm]')
end
colormap(cmap)

    cbh = colorbar;
    if ismember(what,{'phase', 'Phase'}) 
    cbh.Label.String = 'Phase [rad]';
    end
title(titleStr)
ylabel('[mm]')
end