%we are trying to reconstruct the star images
file_path='E:\Chieco\HokeyPokey\pillar data';

pillar_in=imread([file_path '\imagec0_invert.png']);
% imshow(pillar_in)

xyr_locs=dlmread('E:\Chieco\HokeyPokey\pillar data\imagec0_xyr_mat.txt');
r_min=floor(min(xyr_locs(:,3)));
r_max=ceil(max(xyr_locs(:,3)));
box_size=2*r_max+1;

%now we generate the stamps
stamps=zeros(box_size,box_size,r_max-r_min+1);
radius_vec=(r_min:r_max);
%Now we construct our circular mask
% x1 = box_size/2;
% y1 = box_size/2;
%Generate grid with binary mask representing the circle.
% [xx,yy] = meshgrid((1:box_size)-y1,(1:box_size)-x1);
[xx,yy] = meshgrid((-r_max:r_max),(-r_max:r_max));

for i1=1:numel(radius_vec)
    radius=radius_vec(i1);
    mask=zeros(box_size,box_size);
    mask((xx.^2 + yy.^2)<radius^2)=1;
    stamps(:,:,i1)=mask;
%     imwrite(stamps(:,:,i1),['E:\Chieco\HokeyPokey\pillar data\masks\stamp_' num2str(radius) '.png'])
end

%we always want dark backgrounds and bright particles
background_value_max=91;
imi_in=double(pillar_in);
[vertex_locations,vertex_imi_final]=pillar_centers(imi_in,stamps,round(xyr_locs(:,1:2)),background_value_max);

imwrite(vertex_imi_final,[file_path '\mask_reconstruct_from_image.png'])


