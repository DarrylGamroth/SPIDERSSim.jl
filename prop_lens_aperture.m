function [bm, loss] = prop_lens_aperture(bm, fl, dia, snm, varargin)
 
p = inputParser;
p.addParameter('disp',0, @isnumeric);
p.addParameter('power',0.4);
p.addParameter('zoom',4, @isnumeric);
p.parse(varargin{:});

inbm = bm;
bm = prop_lens(bm,fl,snm);
 
 influx = sum(sum(prop_get_amplitude(bm).^2));
 bm   = prop_circular_aperture(bm,dia/2);
 outflux = sum(sum(prop_get_amplitude(bm).^2));
 loss = influx/outflux - 1;

 if p.Results.disp 
     dim = prop_get_gridsize();
     dx = prop_get_sampling(bm);
    x = linspace(-dim/2,dim/2-1,dim)*dx;
    
     intens = prop_get_amplitude(inbm).^2;
     if isnumeric(p.Results.power)
         ima = intens; % .^p.Results.power;
         if p.Results.power ==1
             cmap = cubehelix(512,0.84,0.2,2.66,1);
         else
             cmap = cubehelix(4096,0.84,0.2,2.66,p.Results.power); % more colors because of gamma compression/expansion
         end
     else
         ima = log(intens);
         cmap = cubehelix(512,0.84,0.2,2.66,1);
     end
     
    figure
    imagesc(x*1e3,x*1e3,ima); axis equal tight; zoom(p.Results.zoom)
    colormap(cmap)
    colorbar
    title(sprintf('%s | Dia.= %3.2f mm | Vignetting = %3.2e', snm, dia*1e3, loss))
    xlabel('[mm]')
    ylabel('[mm]')
    % draw lens aperture
    th = 0:pi/100:2*pi;
    xunit = dia/2*1e3 * cos(th);
    yunit = dia/2*1e3 * sin(th);
    hold on
    plot(xunit, yunit, 'r-');
end
 
end