function fft_radial_out=fft_radial_sparse(im1,im2,norm1,norm2)

% Image size and normalization
rows=numel(im1(:,1));
cols=numel(im1(1,:));

% Fourier transform
F1 = fftshift(fft2(fftshift(full(im1))));
clear im1
F2 = fftshift(fft2(fftshift(full(im2))));
clear im2
F = (abs(F1*norm1+F2*norm2).^2);
clear F1
clear F2
% Show spectrum (log)
% imagesc(log(F)); title('Fourier transform (abs log)'); 


freq=(-cols/2:cols/2-1);
[X, Y]=meshgrid(freq,freq); % get coordinates of the power spectrum image
[theta, rho]=cart2pol(X,Y); % equivalent in polar coordinates

%we grid our space in terms of wave vectors
% [ux,uy]=meshgrid(2*pi*(-cols/2:cols/2-1)/cols,2*pi*(-rows/2:rows/2-1)/rows);
% %Convert to radial wave vector
% rcoords = (0:sqrt((2*pi/rows)^2 + (2*pi/cols)^2):sqrt((2*pi*(rows-1)/rows)^2 + (2*pi*(cols-1)/cols)^2));
% % thcoords = linspace(0,2*pi,cols);
% % [ri,thi] = meshgrid(rcoords,thcoords);
% % [x,y] = pol2cart(thi,ri);
% % Fp = interp2(X,Y,F,ux,uy);
% [theta, rho]=cart2pol(ux,uy);


eff=zeros(floor(rows/2)+1,1);

for r=0:rows/2-1
 eff(r+1)=mean(F(find(and(rho>=r,rho<r+1)))); % average power values of all the same freq
end


% freqency spectrum
freq2=0:rows/2; 
% loglog(freq2,eff); title('frequency spectrum','Fontsize',14); axis tight
% xlabel('Freqencies','Fontsize',12); ylabel('Power','Fontsize',12); 
% bark=0;

% % Grid of FFT coordinates
% [rows, cols] = size(F);
% [ux, uy] = meshgrid(([1:cols]-(fix(cols/2)+1))/(cols-mod(cols,2)), ...
%     ([1:rows]-(fix(rows/2)+1))/(rows-mod(rows,2)));
% % % % Convert to polar coordinates
% % % Fr = F; %.* r;
% % imagesc(abs(F)); title('Fourier transform x radius'); 
% % bark=0;
% % rcoords = linspace(min(min(r(r ~= 0))),sqrt(ux(1,1)^2 + uy(1,1)^2),rows);
% rcoords = linspace(0,sqrt(ux(1,1)^2 + uy(1,1)^2),rows);
% thcoords = linspace(0,2*pi,cols);
% [ri,thi] = meshgrid(rcoords,thcoords);
% [x,y] = pol2cart(thi,ri);
% Fp = interp2(ux,uy,F,x,y);
% % imagesc(Fp); title('Fourier transform in polar coordinates'); 
% % bark=0;
% % Sum columns to give 1D projection
% F1D = sum(Fp);
fft_radial_out=[freq2',eff];