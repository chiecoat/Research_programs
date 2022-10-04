function fft_radial_out=fft_radial_large(im1,norm1)

% Image size and normalization
rows=numel(im1(:,1));
cols=numel(im1(1,:));

% Fourier transform
F1 = fftshift(fft2(fftshift(full(im1))));
clear im1
F = (abs(F1*norm1).^2);
% imshow(log(F),[ceil(min(min(log(F)))) ceil(max(max(log(F))))]); colormap(jet); colorbar
% clear F1

% Show spectrum (log)
% imagesc(log(F)); title('Fourier transform (abs log)'); 


freq=(-cols/2:cols/2-1);
[X, Y]=meshgrid(freq,freq); % get coordinates of the power spectrum image
[theta, rho]=cart2pol(X,Y); % equivalent in polar coordinates

% eff=zeros(floor(rows/2)+1,1);
% 
% % for r=0:rows/2-1
% for r=0:11
%     vals=F(find(and(rho>=r,rho<r+1)));
%     eff(r+1)=mean(F(find(and(rho>=r,rho<r+1)))); % average power values of all the same freq
%     if r>0
%         semilogy(zeros(numel(vals),1)+r,vals,'o');
%         hold on
%     end
% end
% plot((1:11),eff(2:12),'-+','LineWidth',4)


r_vec=unique(rho);
eff=zeros(numel(r_vec),1);

for r=1:numel(r_vec)
    vals=F(find(rho==r_vec(r)));
    eff(r)=mean(F(find(rho==r_vec(r)))); % average power values of all the same freq
%     if r>0
%         semilogy(zeros(numel(vals),1)+r_vec(r),vals,'o');
%         hold on
%     end
end
% plot(r_vec(2:41),eff(2:41),'-+','LineWidth',4)

freq2=0:rows/2; 

fft_radial_out=[freq2',eff];