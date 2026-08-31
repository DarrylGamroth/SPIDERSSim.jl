function out = cdiOTF(bm,d,sep, pa, varargin)
% d : pupil dia in pupil dia (1)
% sep : side lobe position in pupil ddia(1.6)
% pa : position angle of the side lobe

p = inputParser;
p.addParameter('disp',0, @isnumeric);
p.parse(varargin{:});

np = prop_get_gridsize();
dx = prop_get_sampling(bm);
dxr = prop_get_sampling_radians(bm);
lamb = prop_get_wavelength(bm);

phd = d/0.9*1/43*24/14; % pinhole ddia in cropped pup dia

dxp = lamb/bm.diam/dxr/np; % sampling of the OTF plan in pupil diameter

xc = sep*cos(pa)/dxp;
yc = sep*sin(pa)/dxp;
maskPh0 = prop_ellipse(bm,phd/2/dxp*dx, phd/2/dxp*dx); % mask for pinhole
maskPh0b = prop_ellipse(bm,2*phd/dxp*dx, 2*phd/dxp*dx) - maskPh0; % mask for pinhole surrondings

mask0 = prop_ellipse(bm,d/dxp*dx, d/dxp*dx); % mask for central lobe
mask1 = prop_ellipse(bm,d/2/dxp*dx, d/2/dxp*dx,'xc',xc*dx,'yc',yc*dx); % mask for sidelobe
ima = (prop_get_amplitude(bm)).^2;
otf = fftshift(fft2(ima));

p0 = abs(ifft2(otf.*mask0));
p1c = circshift(otf.*mask1,round([-yc -xc]));
p1 = abs(ifft2(p1c));
vis2D = 0*p0;
vis2D(p0>0) = 2*p1(p0>0)./p0(p0>0);

% compute Ir
tmp = maskPh0b.*otf;
phBg = median(tmp(:)); % median level around pinhole MTF
ir = abs(ifft2((otf-1*phBg).*maskPh0));
out.otf = otf;
out.mtf = abs(otf);
out.ir = ir;
out.iminus = p1;
out.i0 = p0;
out.ip = p0-p1.^2./ir-ir;
out.ip(out.ip<0) = 0;

out.ip2 = 0*p0;
out.ip2(vis2D>0) = p0(vis2D>0).*(1-vis2D(vis2D>0));
out.ip2(out.ip2<0) = 0;
% out.ip = p0-p1;
out.vis2D = vis2D;
out.vis = 2*sum(p1(:))/sum(p0(:));
out.extinction1 = out.i0./out.ip;
out.extinction2 = out.i0./out.ip2;

if p.Results.disp
    % np = prop_get_gridsize();
    x = linspace(-np/2,np/2-1,np);
    % fx = 1/dx * (-np/2:np/2-1)/np;
    fx = x.*lamb/bm.diam/dxr/np;
    % [fxx,fyy]=meshgrid(fx,fx);
  
    figure
    subplot(3,3,1)
    imagesc(fx,fx,out.mtf.^0.3); axis equal tight; zoom(2)
    title(sprintf('SCC MTF | Fringe Vis. = %0.3f', out.vis))
    xlabel('[pupil dia.]')
    ylabel('[pupil dia.]')
    
    subplot(3,3,2)
    imagesc(x*dx*1e3,x*dx*1e3,out.ir.^0.3); axis equal tight; zoom(1)
    title('Pinhole PSF')
    xlabel('[mm]')
    ylabel('[mm]')
     
    subplot(3,3,3)
    imagesc(x*dx*1e3,x*dx*1e3,out.i0.^0.3); axis equal tight; zoom(1)
    title('I0')
    xlabel('[mm]')
    ylabel('[mm]')
    
    subplot(3,3,4)
    imagesc(x*dx*1e3,x*dx*1e3,out.iminus.^0.3); axis equal tight; zoom(1)
    title('I-')
    xlabel('[mm]')
    ylabel('[mm]')
    
    subplot(3,3,5)
    imagesc(x*dx*1e3,x*dx*1e3,out.vis2D,[0 1]); axis equal tight; zoom(1)
    title(sprintf('SCC Fringe Visibility Map | Average Fringe Vis. = %0.3f', out.vis))
    xlabel('[mm]')
    ylabel('[mm]')
    
    subplot(3,3,6)
    imagesc(x*dx*1e3,x*dx*1e3,out.ip.^0.3); axis equal tight; zoom(1)
    title('I planet regular')
    xlabel('[mm]')
    ylabel('[mm]')
     
    subplot(3,3,7)
    imagesc(x*dx*1e3,x*dx*1e3,out.ip2.^0.3); axis equal tight; zoom(1)
    title('I planet using Vis')
    xlabel('[mm]')
    ylabel('[mm]')   
    
    subplot(3,3,8)
    imagesc(x*dx*1e3,x*dx*1e3,out.extinction1); axis equal tight; zoom(1)
    title('CDI extinction map (regular)')
    xlabel('[mm]')
    ylabel('[mm]')   
        
    subplot(3,3,9)
    imagesc(x*dx*1e3,x*dx*1e3,out.extinction2); axis equal tight; zoom(1)
    title('CDI extinction map (using VIS)')
    xlabel('[mm]')
    ylabel('[mm]')   
%     
%     
% figure
%     imagesc(x*dx*1e3,x*dx*1e3,1./(1-out.vis2D), [0 100]); axis equal tight; zoom(1)
%     title('CDI extinction map : 1/(1-VIS)')
%     xlabel('[mm]')
%     ylabel('[mm]') 
    
      figure
subplot(1,3,1)
    imagesc(fx,fx,out.mtf.^0.3.*(mask0+mask1)); axis equal tight; zoom(2)
    title(sprintf('SCC MTF x Masks| Fringe Vis. = %0.3f', out.vis))
    xlabel('[cycle/m]')
    ylabel('[cycle/m]')
    subplot(1,3,2)
    imagesc(fx,fx,out.mtf.^0.3.*(maskPh0)); axis equal tight; zoom(2)
    title(sprintf('SCC MTF x Pinhole Mask| Fringe Vis. = %0.3f', out.vis))
    xlabel('[cycle/m]')
    ylabel('[cycle/m]')
        subplot(1,3,3)
    imagesc(fx,fx,out.mtf.^0.3.*(maskPh0b)); axis equal tight; zoom(2)
    title(sprintf('SCC MTF x Pinhole surrounding Mask| Fringe Vis. = %0.3f', out.vis))
    xlabel('[cycle/m]')
    ylabel('[cycle/m]')
end


end