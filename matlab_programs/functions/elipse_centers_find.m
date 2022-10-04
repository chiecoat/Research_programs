%we are trying to reconstruct the star images
file_path='E:\Chieco\HokeyPokey\ellipse data';

x_offset=500;
y_offset=180;

my_try_bin=imread([file_path '\erode_binary_shell_img0_cropFrom_x' num2str(x_offset) '_y' num2str(y_offset) '.png']);
my_try_bin=double(my_try_bin);
imshow(my_try_bin)

watershed_mat=watershed(my_try_bin,4);
rgb = label2rgb(watershed_mat,'jet',[.5 .5 .5]);
imshow(rgb);
hold on

imwrite(rgb,'E:\Chieco\HokeyPokey\ellipse data\watershed_ellipses.png')

n_rows=numel(watershed_mat(:,1));
n_cols=numel(watershed_mat(1,:));
imi_pad=uint16(zeros(n_rows+77*2,n_cols+77*2)+1);
imi_pad(77:77+n_rows-1,77:77+n_cols-1)=watershed_mat;
imwrite(imi_pad,[file_path '\erode_ids_img0_cropFrom_x' num2str(x_offset) '_y' num2str(y_offset) '.png']);
writematrix(imi_pad,[file_path '\erode_ids_img0_cropFrom_x' num2str(x_offset) '_y' num2str(y_offset) '.csv'], 'delimiter',',')

%now we want to get data from the watershed images
try_region=regionprops(imi_pad,'Area','Centroid','Circularity','Eccentricity','Perimeter');
areas=[try_region.Area];
n_regions=numel(areas);
circ=[try_region.Circularity];
ecc=[try_region.Eccentricity];
perim=[try_region.Perimeter];
%we need to identify the particle centers
cens=[try_region.Centroid];
x_cens=cens((1:2:2*n_regions-1));
y_cens=cens((2:2:2*n_regions));
%we put all the data tgether in to one matrix
bub_mat=[x_cens(2:end)',y_cens(2:end)',areas(2:end)',...
    circ(2:end)',ecc(2:end)',perim(2:end)',(2:n_regions)'];
bub_mat=bub_mat(bub_mat(:,3)<1E6,:);

%the first point is the interstitial
bub_cens_keep=bub_mat(bub_mat(:,3)>1E2,:);
figure(1); 
plot(bub_cens_keep(:,1),bub_cens_keep(:,2),'ok','Linewidth',2)
figure(2)
semilogy(bub_mat(:,4),bub_mat(:,3),'or')
figure(3)
semilogx(bub_mat(:,3),(bub_mat(:,6).^2)./bub_mat(:,3),'or')
close Figure 1 
close Figure 2 
close Figure 3 

im_reg=imread([file_path '\raw data\img00000.tif']);
imshow(imi_pad/255);
hold on
% plot(bub_cens_keep(:,1)+x_offset,bub_cens_keep(:,2)+y_offset,'or','Linewidth',2)
plot(bub_cens_keep(:,1),bub_cens_keep(:,2),'or','Linewidth',2)

file_write=[file_path '\ellipse_centers_img0_crop.txt'];
dlmwrite(file_write,[bub_cens_keep(:,1),bub_cens_keep(:,2)],'delimiter',',')
keyboard