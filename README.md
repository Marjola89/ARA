# ARA

Mutual Information was calculated from [amis corrs] = ami(xy,nBins,nLags) by Durga Lal Shrestha.
https://uk.mathworks.com/matlabcentral/fileexchange/7936-ami-and-correlation?focused=6141194&tab=function

Call from the main programm and plot as:

N=reshape(Est',nbins,nbins);
pts = linspace(min(v), max(v), nbins);

% Create Gaussian filter matrix:
[xG, yG] = meshgrid(-(nbins-0.5)/2:(nbins-0.5)/2);
sigma = 2.5;
g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
g = g./sum(g(:));

figure;
imagesc(pts, pts, N);xlabel('v');ylabel('w');
set(gca,'FontSize',15);
caxis(gca, 'manual');
colormap(jet);
colorbar;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
axis tight;
set(gca,'Color',[0.04 0.04 0.52]);

%
