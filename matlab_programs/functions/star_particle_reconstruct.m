%we are trying to reconstruct the star images
file_path='E:\Chieco\HokeyPokey';

my_try_bin=imread([file_path '\stars_im_erode1_edges.png']);
imshow(my_try_bin)

watershed_mat=watershed(my_try_bin,8);
rgb = label2rgb(watershed_mat,'jet',[.5 .5 .5]);
imshow(rgb);
hold on



%now we want to get data from the watershed images
try_region=regionprops(watershed_mat,'Area','Centroid','Circularity','Eccentricity','Perimeter');
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
bub_cens_keep=bub_mat(bub_mat(:,3)>50,1:2);
figure(1); 
plot(bub_cens_keep(:,1),bub_cens_keep(:,2),'ok','Linewidth',2)
figure(2)
semilogy(bub_mat(:,4),bub_mat(:,3),'or')
figure(3)
semilogx(bub_mat(:,3),(bub_mat(:,6).^2)./bub_mat(:,3),'or')
close Figure 1 
close Figure 2 
close Figure 3 

%Now we read in the structuring element
struct=imread([file_path '\star_structure_element_0.png']);
%this is for visualizing the structuring element
watershed_struct=watershed(struct,8);
rgb = label2rgb(watershed_struct,'jet',[.5 .5 .5]);
imshow(rgb);
watershed_struct(watershed_struct==0)=2;
watershed_struct(watershed_struct==1)=0;
watershed_struct(watershed_struct>1)=1;
struct_fin=watershed_struct;
%this is a weirf quirk that we will fix later. need the box to have odd
%number of rows and columns
struct_fin=[struct_fin,zeros(52,1)]; 
struct_fin=[struct_fin;zeros(1,53)];
box_size=53;

n_rows=numel(struct_fin(:,1));
n_cols=numel(struct_fin(1,:));
%this is known from the symmetry of the partilcle
n_thetas=72;
s_rot=zeros(n_rows,n_cols,n_thetas);
%this would be for a generic structuring element angle
% t_start=-vertex_loc_clean(3); t_fin=120-vertex_loc_clean(3)-1;

%right now we have a strusturing element that starts at theta=0
t_start=0; t_fin=71;
for theta=t_start:t_fin
    s_new=imrotate(struct_fin,theta);
    sx_cen=floor(numel(s_new(1,:))/2)+1; sy_cen=floor(numel(s_new(:,1))/2)+1;
    s_new=s_new(sy_cen-floor(box_size/2):sy_cen+floor(box_size/2),...
        sx_cen-floor(box_size/2):sx_cen+floor(box_size/2));
    s_rot(:,:,theta+1)=s_new;
%     imwrite(s_rot(:,:,theta+1),[file_path '\structure elements angles\star_structure_element_' num2str(theta) '.png'])
end


background_value=0;
imi_in=imread([file_path '\stars_im_machine.png']);%double(watershed_mat);
imi_in=double(imi_in);
[vertex_locations,vertex_imi_final]=star_watershed_struct_2d(imi_in,s_rot,bub_cens_keep,background_value);

imwrite(vertex_imi_final,[file_path '\star_reconstruct_final_from_ML.png'])


