function out = prop_psfContrast(bm, ref)

psf = (prop_get_amplitude(bm)).^2;
np = prop_get_gridsize();
dx = prop_get_sampling(bm);
dxr = prop_get_sampling_radians(bm);
lamb = prop_get_wavelength(bm);
x = linspace(-np/2,np/2-1,np);
[xx yy] = meshgrid(x,x);
rho = sqrt(xx.^2+yy.^2); % radius in pixels

r = 0:np/2-1; % radial space (x>=0)
resr = r*dxr/lamb*bm.diam; % axis in resel

quan = [0.25 0.5 0.75]; % percent quantile to compute
cont = zeros(numel(quan), numel(r));

for r = 0:np/2-1 % in pixels
    cont(:,r+1) = quantile(psf(abs(rho-r)<=0.5), quan);

end

if ref % normalize contrat curve to max of non-coro PSF
    refPsf = importdata('nonCoronPSF.mat');
    cont = cont./max(refPsf.psf2D(:));
    refCont = refPsf.cont./max(refPsf.psf2D(:));
end
figure
jbfill(resr,cont(3,:),cont(1,:),'r','w',1,0.2); % shade area between 25% and 75% curves
hold on
plot(resr,cont(2,:),'r-','lineWidth', 2);
if ref
    jbfill(resr,refCont(3,:),refCont(1,:),'b','w',1,0.2); % shade area between 25% and 75% curves
    hold on 
    plot(refPsf.resel,refCont(2,:),'b-', 'lineWidth',1);
end
set(gca,'Yscale', 'log'); % jbfill undo the log scale
% xlim([0 25]);
grid
xlabel('Radius [\lambda/D]');
ylabel('Contrast');
legend('Coro PSF 25% to 75%-tile contrast', 'Coro PSF median contrast', ...
    'Non- Coro PSF 25% to 75%-tile contrast','Non-Coro PSF median contrast');

out.resel = resr;
out.percentile = quan;
out.cont = cont;
end