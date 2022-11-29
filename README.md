# ARA

The main aim of this analysis is to perform and assess the variability of a quasi-periodic signal as a form of an attractor. Provides a visual representation of the multidimensional dynamics of a time series from a two-dimensional attractor with features like density and symmetry. This method is used for the analysis of microvascular blood flow signals and details can be found in the conference paper below: 

M. Thanaj, A. J. Chipperfield, and G. F. Clough, “Attractor Reconstruction Analysis for Blood Flow Signals,” Conference: 41st International Engineering in Medicine and Biology Conference (EMBC 2019).





Mutual Information was calculated from [amis corrs] = ami(xy,nBins,nLags) by Durga Lal Shrestha.
https://uk.mathworks.com/matlabcentral/fileexchange/7936-ami-and-correlation?focused=6141194&tab=function

Call from the main programm as:
```
[Est,Dt,mean_Angles,Angle_rotation,Angular_spread,Radius,Period,time_lag,u,v,w] = uvw_analysis(x_data, Fs, nbins);
maximum_Density=max(Dt);

and plot as:

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

```
