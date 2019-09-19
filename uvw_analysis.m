 function [Est,Dt,angp,thetar,thetas,radius,T,tau,u,v,w] = uvw_analysis(x_w, Fs, nbins)
% x_w    :  The window of the data of observation
% Fs     :  Sampling frequency
% nbins  :  Number of bins for estimating the kernel density

Mlag=floor(length(x_w)/2);
[amis, ~, ~] = ami(x_w,20,Mlag);estLag = find(diff(amis)>0,1);% compute the time lag
[pksh,lcsh] = findpeaks(amis);short = mean(diff(lcsh))/(Fs);
if ~(isnan(short))
    [~,lclg] = findpeaks(amis,'MinPeakHeight',max(amis(estLag:end))/4,'MinPeakDistance',ceil(short)*(Fs));
    tao = find(lclg > floor(Fs/2) & lclg < floor(Fs*1.6),1);
    
    if isempty(tao)
        T=estLag*2;
    else
        T=lclg(tao);
    end
else
    T=floor(3*estLag);
end
tau=floor(T/3);
[u,v,w]=uvw_rec(x_w,tau);
[v_rot, w_rot] = rotate_plane(v, w,2);
[v_rot4, w_rot4] = rotate_plane(v, w,4);
% maximum ksdensity 
[Est, ~] = densimage_calc(v, w,nbins);
[Estrot, ~] = densimage_calc(v_rot', w_rot',nbins);
[Estrot4, ~] = densimage_calc(v_rot4', w_rot4',nbins);
% based on the paper
Dt=(Est+Estrot+Estrot4)/3; dt=max(Dt); dtmin=min(Dt);
dest=max(Est);

[theta,thetar,radius]=angle_calc(v,w,x_w,tau); 
[theta_rot,~,~]=angle_calc(v_rot,w_rot,x_w,tau);
[theta_rot4,~,~]=angle_calc(v_rot4,w_rot4,x_w,tau);
thet=[theta theta_rot theta_rot4];
angp=mean(thet);thetas=abs(mean(diff(thet)));

end

% rotate u, v and w
function [v_rotated, w_rotated] = rotate_plane(v, w, th)
% rotation th = 2*pi/3 or th = 4*pi/3
p = [v w]';
% choose a point which will be the center of rotation
kerSize = 6 ;
v_smz = conv( v, ones(1, kerSize)/kerSize, 'same' ) ;sv = sign( v_smz ) ;
tIdv = find( sv(1:end-1) .* sv(2:end) < 1 ) ;
w_smz = conv( w, ones(1, kerSize)/kerSize, 'same' ) ;sw = sign( w_smz ) ;
tIdw = find( sw(1:end-1) .* sw(2:end) < 1 ) ;
x_center = v(tIdv(2));
y_center = w(tIdw(2));
% create a matrix which will be used later in calculations
center = repmat([x_center; y_center], 1, length(v));
% define a 60 degree counter-clockwise rotation matrix
theta = (th*pi)/3;       % pi/3 radians = 60 degrees
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% do the rotation...
s = p - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
% this can be done in one line as:
% vo = R*(v - center) + center
% pick out the vectors of rotated x- and y-data
v_rotated = vo(1,:);
w_rotated = vo(2,:);
% make a plot
% figure;
% plot(v, w, 'k-', v_rotated, w_rotated, 'r-', x_center, y_center, 'bo');
% axis tight

end

% Kernel density estimation
function [Est, XY] = densimage_calc(v, w,nbins)
% nbins=47;
% fact=factor(length(v));mgrid1=fact(end);mgrid2=(length(v))/fact(end);
gridx1 = linspace(min(v)*1.5,max(v)*1.5,nbins);
gridx2 = linspace(min(w)*1.5,max(w)*1.5,nbins);
[x1,x2] = meshgrid(gridx1, gridx2);
% [x1,x2] = meshgrid(-(nbins-0.5)/2:(nbins-0.5)/2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
% hopt=1;
[Esto,XYo,bw] = ksdensity([v w],xi);%hopt=bw/2.8;

% optimase the kernel bandwidth
htry=(0.3:0.05:1.5).*bw(1); htry2=(0.3:0.05:1.5).*bw(2);
for i=1:length(htry)
    moh(i)=lscvscore(v,htry(i));
    moh2(i)=lscvscore(w,htry2(i));
end
hopt1=htry(find(moh <=min(moh)));hopt2=htry2(find(moh2 <=min(moh2)));
 hopt=mean([hopt1 hopt2])/1.5; % optimised bandwidth
[Est,XY]=ksdensity([v w],xi,'Bandwidth',hopt);
 
end
% optimise the bandwidth of the kernel smoothing window using
% the least-square cross validation approach with Gaussian kernels. 
function moh=lscvscore(x,h)
% x is the input data
% h is the bandwidth
n=length(x);
sqterm=0;
xterm=0;
for i=1:n
    sqterm=sqterm+sum(kg((x-x(i))/(sqrt(2)*h)))/sqrt(2);   %  The sqrt(2) factors are for 
    % the K(2) term which is
    % equivalent to a kernel with variance of 2 that is the convolution of 2
    % gaussian kernels.  
    xterm=xterm+sum(kg((x-x(i))/h));
end
sqterm=sqterm/(n*n*h);
xterm=2*(xterm/(n*(n-1))-kg(0)/(n-1))/(h);
moh=sqterm-xterm;
end
% Gaussian Kernel 
function k=kg(x)
k=exp(-x .* x ./2)/sqrt(2*pi);
end
% find the angle of two sides of the triangle
function [theta,thetar,r]=angle_calc(v,w,x_w,la3)

% estimate angle of rotation and angular spread
% [hv,pv]=findpeaks(v,'MinPeakHeight',max(v)/8,'MinPeakDistance',2*T3/3);
[hvn,~]=findpeaks(-v,'MinPeakDistance',2*la3);% amplitude of negative v
[hw,~]=findpeaks(w,'MinPeakDistance',2*la3);% amplitude of positive w
[hwn,~]=findpeaks(-w,'MinPeakDistance',2*la3);% amplitude of negative w
[hx,~]=findpeaks(x_w,'MinPeakDistance',2*la3,'MinPeakHeight',mean(x_w));% amplitude of original signal
hamp=hx-mean(x_w); 
xdetr=detrend(x_w); % detrend to find negative peaks 
[hnx,pnx]=findpeaks(-xdetr,'MinPeakDistance',2*la3,'MinPeakHeight',mean(xdetr));
hny=x_w(pnx);hamp_n=mean(hx)-mean(hny);
% hy=peak2peak(x_w); % wrong because it is not a stable signal
rwn=mean(abs(hwn));rw=mean(abs(hw));
rvn=mean(abs(hvn));%rv=mean(abs(hv));
% the equilateral triangle centred on the origin with the bottom edge given as 
% w = - h / (2*sqrt(2)) then the angle of rotation is defined by w = - h*sec(theta) / (2*sqrt(2)) 
% therefore, theta = asec(- w*2*sqrt(2) / h).
% hwnew = hwn(1:30);
w_new=-(w*2*sqrt(2));[hw_new,~]=findpeaks(w_new,'MinPeakDistance',2*la3);% amplitude of positive w new
X = rdivide(mean(hw_new),(hamp_n)); % calculate the median values of the amplitudes 
thetar=abs((asec(X)));%measuring the rotation
r=max((sqrt(3)*hamp_n)/(2*sqrt(2))); % radius of the circle centred at (0,0) encloses 95% of A
theta=asec((rwn+rw)/rvn);%measuring the angle of triangles
end

% compute u, v and w
function [u,v,w]=uvw_rec(x_w,lag3)
em=3;
N1=length(x_w);M1=N1-(em-1)*lag3(1,1);Y3=zeros(M1,em);

for i=1:em
    Y3(:,i)=x_w((1:M1)+(i-1)*lag3)';
end
xt=Y3(:,3);yt=Y3(:,2);zt=Y3(:,1);
% remove the baseline vartiation by projecting it in a new plane (v,w)
% prependicular to (1,1,1) vector.
u=(1/sqrt(3))*(xt+yt+zt);
v=(1/sqrt(6))*(xt+yt-(2*zt));
w=(1/sqrt(2))*(xt-yt);
end
