readPath='E:\Chieco\Hyperuniformity\HU foams\track data';
im_read=[125,250,500];

%we read in the data and show it collapses for A/<A>
figure(1)
for i1=1:numel(im_read)
    bub_lengths_mat=[-1,-1,-1];
    im_num_str=num2str(im_read(i1));
    vert_mat=dlmread([readPath '\vertices\vertices imi_' im_num_str '_final.txt']);
    bub_mat=dlmread([readPath '\bubbles\bubbles imi_' im_num_str '_final.txt']);    
    stats_mat=dlmread([readPath '\stats\stats_bub_nxy_pace imi_' im_num_str '_final.txt']);
    %we need to remove the bubbles that connect to the edge. We do this by
    %identifying vertices that touch cell=1 and remove any bubbles with
    %that vertex from the list
    cell_ids=vert_mat(:,end-2:end);
    id_list=find(cell_ids==1);
    [vert_rows,vert_cols]=ind2sub(size(cell_ids),id_list);
    bub_ids=bub_mat(:,1);
    vert_ids_mat=bub_mat(:,5:end);
    sz_of_verts=size(vert_ids_mat);
    for i2=1:numel(vert_rows)
        bub_id_list=find(vert_ids_mat==vert_rows(i2));
        [bub_rows,bub_cols]=ind2sub(sz_of_verts,bub_id_list);
        bub_ids(bub_rows)=-1;        
    end
    bub_keep=stats_mat(bub_ids>0,:);
    if i1==1
        bub_keep_tot=bub_keep;
    else
        bub_keep_tot=[bub_keep_tot;bub_keep];
    end
    norm=numel(bub_keep(:,1));
    norm_area=sum(bub_keep(:,6));
    %we find p_of_n by looping through the number of sides
    n_sides=unique(bub_keep(:,1));
    sides_tot=zeros(numel(n_sides),5);
    for i3=1:numel(n_sides)        
        bub_sides=bub_keep(bub_keep(:,1)==n_sides(i3),:);
        pOfn=numel(bub_sides(:,1))/norm;
        del_pOfn=pOfn/sqrt(norm);
        FOfn=sum(bub_sides(:,6))/norm_area;
        del_FOfn=std(bub_sides(:,6)/norm_area);
        sides_tot(i3,:)=[n_sides(i3),pOfn,del_pOfn,FOfn,del_FOfn];        
    end
    figure(1)
    histogram(bub_keep(:,1),'Normalization','pdf');
    plot(sides_tot(:,1),sides_tot(:,2),'o')
    hold on
    figure(2)
    semilogy(sides_tot(:,1),sides_tot(:,4),'LineWidth',2);
    hold on   
    cdf_imi_name=['E:\Chieco\Hyperuniformity\HU foams\distribution data\foam data\foam_pOfn_FOfn imi_' im_num_str '.txt'];
    dlmwrite(cdf_imi_name,sides_tot,'Newline','pc')
end

norm_tot=numel(bub_keep_tot(:,1));
norm_area_tot=sum(bub_keep_tot(:,6));
%we find p_of_n by looping through the number of sides
n_sides_tot=unique(bub_keep_tot(:,1));
sides_scale=zeros(numel(n_sides_tot),5);
for i3=1:numel(n_sides_tot)
    bub_sides_tot=bub_keep_tot(bub_keep_tot(:,1)==n_sides_tot(i3),:);
    pOfn=numel(bub_sides_tot(:,1))/norm_tot;
    del_pOfn=pOfn/sqrt(norm_tot);
    FOfn=sum(bub_sides_tot(:,6))/norm_area_tot;
    del_FOfn=std(bub_sides_tot(:,6)/norm_area_tot)/sqrt(numel(im_read));
    sides_scale(i3,:)=[n_sides_tot(i3),pOfn,del_pOfn,FOfn,del_FOfn];
end
figure(1)
histogram(bub_keep_tot(:,1),'Normalization','pdf');
plot(sides_tot(:,1),sides_tot(:,2),'o')
figure(2)
semilogy(sides_scale(:,1),sides_scale(:,4),'LineWidth',2);

total_name='E:\Chieco\Hyperuniformity\HU foams\distribution data\foam data\foam_pOfn_FOfn total.txt';
dlmwrite(total_name,sides_scale,'Newline','pc')


% %new we perform the same analysis but for the voronoi patterns.
% readPath='C:\Users\Seyyed\Desktop\Files ATC\hyperuniformity foams\particle positions\voronoi\matrix to convert';
% pattern={'Halton','Einstein_square_delta026','Einstein_square','poisson'};
% 
% for i1=1:numel(pattern)
%     centroid_dist_mat=[-1,-1,-1];
%     points_dist_mat=[-1,-1,-1];
%     pat=pattern{i1};
%     cell_mat=dlmread([readPath '\' pat ' ID_NA_xyp_xyc_vertID N_cells_500000 run_1.txt']);
%     vert_mat=dlmread([readPath '\' pat ' vertex_xyID N_cells_500000 run_1.txt']);
%     %we need to remove the bubbles that connect to the edge. We find the
%     %edge vertices and find which bubbles they are on
%     vert_left=vert_mat(vert_mat(:,1)==0,3);
%     vert_right=vert_mat(vert_mat(:,1)==1,3);
%     vert_bot=vert_mat(vert_mat(:,2)==0,3);
%     vert_top=vert_mat(vert_mat(:,2)==1,3);
%     bad_verts=unique([vert_left;vert_right;vert_top;vert_bot]);
%     %now we want to find the bubbles that have a bad vertex
%     bub_ids=cell_mat(:,1);
%     vert_ids_mat=cell_mat(:,8:end);
%     sz_of_verts=size(vert_ids_mat);
%     numel(bad_verts)
%     for i2=1:numel(bad_verts)
%         bub_id_list=find(vert_ids_mat==bad_verts(i2));
%         [bub_rows,bub_cols]=ind2sub(sz_of_verts,bub_id_list);
%         bub_ids(bub_rows)=-1;
%     end
%     bub_keep=cell_mat(bub_ids>0,:);
%     norm=numel(bub_keep);
%     %we find p_of_n by looping through the number of sides
%     n_sides=unique(bub_keep(:,2));
%     sides_tot=zeros(numel(n_sides),3);
%     for i3=1:numel(n_sides)        
%         bub_sides=bub_keep(bub_keep(:,2)==n_sides(i3),:);
%         sides_tot(i3,:)=[n_sides(i3),numel(bub_sides(:,2))/norm,sum(bub_sides(:,3))];        
%     end
%     figure(1)
%     histogram(bub_keep(:,2),'Normalization','pdf');
%     hold on
%     figure(2)
%     semilogy(sides_tot(:,1),sides_tot(:,3),'LineWidth',2);
%     hold on    
% 
%     writePath='C:\Users\Seyyed\Desktop\Files ATC\hyperuniformity foams\distributions';
%     size(cdf_write)
%     dlmwrite([writePath '\' pat ' pOfn_Fofn N_cells_500000.txt'],cdf_write)
% end