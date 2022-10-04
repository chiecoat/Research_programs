readPath='F:\Chieco\Hyperuniformity\HU foams\track data';
im_read=[125,250,500];



for i1=1:numel(im_read)
    bub_lengths_mat=[-1,-1,-1];
    im_num_str=num2str(im_read(i1));
    bub_mat=dlmread([readPath '\bubbles\bubbles imi_' im_num_str '_final.txt']);
    n_bubs=numel(bub_mat(:,1));
    vert_mat=bub_mat(:,5:end);
    v_sz=size(vert_mat);
    for i2=1:n_bubs
        id=bub_mat(i2,1);
        verts=vert_mat(i2,vert_mat(i2,:)>0);
        n_verts=numel(verts);
        bub_id_mat=zeros(n_verts,3)-1;
        for i3=1:n_verts
            vert_index=find(vert_mat==verts(i3));
            [row,col]=ind2sub(v_sz,vert_index);
            bub_id_mat(i3,1:numel(row))=row';
        end
        bub_id_vec=reshape(bub_id_mat,[1,3*n_verts]);
        bub_tot_ids=unique(bub_id_vec(bub_id_vec>0));
        bub_neighbors=bub_tot_ids(bub_tot_ids~=id);
        for i4=1:numel(bub_neighbors)
            id_pair=sort([id,bub_neighbors(i4)]);
            if numel(bub_lengths_mat(and(bub_lengths_mat(:,1)==id_pair(1),bub_lengths_mat(:,2)==id_pair(2)),:))>1
               continue
            end
            dist=sqrt((bub_mat(id_pair(1),3)-bub_mat(id_pair(2),3))^2+...
                       (bub_mat(id_pair(1),4)-bub_mat(id_pair(2),4))^2);
            bub_lengths_mat=[bub_lengths_mat;id_pair,dist]; 
        end
    end
    bub_lengths_mat=bub_lengths_mat(2:end,:);
    figure(1)
    hold on;
    histogram(bub_lengths_mat(:,3)/mean(bub_lengths_mat(:,3)));
    if i1==1
        bub_dist_vec=bub_lengths_mat(:,3)/mean(bub_lengths_mat(:,3));
    else
        bub_dist_vec=[bub_dist_vec;bub_lengths_mat(:,3)/mean(bub_lengths_mat(:,3))];
    end    
    cdf_imi_name=['F:\Chieco\Hyperuniformity\HU foams\distribution data\foam data\foam_bub_distance_CDF imi_' im_num_str '.txt'];
    cdf_imi=[0,0;sort(bub_lengths_mat(:,3)/mean(bub_lengths_mat(:,3))),(1:numel(bub_lengths_mat(:,3)))'];
    dlmwrite(cdf_imi_name,[cdf_imi(:,1),1-cdf_imi(:,2)/numel(cdf_imi(:,2))],'Newline','pc')
end

cdf_bub_dist=[0,0;sort(bub_dist_vec),(1:numel(bub_dist_vec))'];


figure(2)
semilogy(cdf_bub_dist(:,1),1-cdf_bub_dist(:,2)/numel(cdf_bub_dist(:,2)),'b','LineWidth',2);
hold on
xs=(0:0.01:6);
ys=exp(-(gamma(1+1/alpha(1))*xs).^alpha(1));
plot(xs,ys,'--r','LineWidth',2)
plot(xs,exp(-xs),':k','LineWidth',2)

total_name='F:\Chieco\Hyperuniformity\HU foams\distribution data\foam data\foam_bub_distance_CDF total.txt';
dlmwrite(total_name,[cdf_bub_dist(:,1),1-cdf_bub_dist(:,2)/numel(cdf_bub_dist(:,2))],'Newline','pc')


% readPath='F:\Chieco\Hyperuniformity\HU foams\voronoi\';
% pattern={'Halton','Einstein_square_delta025','Einstein_square','poisson'};
% 
% for i1=1:numel(pattern)
%     centroid_dist_mat=[-1,-1,-1];
%     points_dist_mat=[-1,-1,-1];
%     pat=pattern{i1};
%     cell_mat=dlmread([readPath '\voro_centroids\all info list\' pat ' ID_NA_xyp_xyc_vertID N_cells_5000 run_1 .txt']);
%     n_bubs=numel(cell_mat(:,1));
%     vert_mat=cell_mat(:,8:end);
%     v_sz=size(vert_mat);
%     for i2=1:n_bubs
%         id=cell_mat(i2,1);
%         verts=vert_mat(i2,vert_mat(i2,:)>0);
%         n_verts=numel(verts);
%         bub_id_mat=zeros(n_verts,3)-1;
%         for i3=1:n_verts
%             vert_index=find(vert_mat==verts(i3));
%             [row,col]=ind2sub(v_sz,vert_index);
%             bub_id_mat(i3,1:numel(row))=row';
%         end
%         bub_id_vec=reshape(bub_id_mat,[1,3*n_verts]);
%         bub_tot_ids=unique(bub_id_vec(bub_id_vec>0));
%         bub_neighbors=bub_tot_ids(bub_tot_ids~=id);
%         for i4=1:numel(bub_neighbors)
%             id_pair=sort([id,bub_neighbors(i4)]);
%             if numel(centroid_dist_mat(and(centroid_dist_mat(:,1)==id_pair(1),centroid_dist_mat(:,2)==id_pair(2)),:))>1
%                continue
%             end
%             dist_cens=sqrt((cell_mat(id_pair(1),4)-cell_mat(id_pair(2),4))^2+...
%                            (cell_mat(id_pair(1),5)-cell_mat(id_pair(2),5))^2);
%             dist_points=sqrt((cell_mat(id_pair(1),6)-cell_mat(id_pair(2),6))^2+...
%                            (cell_mat(id_pair(1),7)-cell_mat(id_pair(2),7))^2);
%             centroid_dist_mat=[centroid_dist_mat;id_pair,dist_cens];
%             points_dist_mat=[points_dist_mat;id_pair,dist_points];
%         end
%     end
%     centroid_dist_mat=centroid_dist_mat(2:end,:);
%     points_dist_mat=points_dist_mat(2:end,:);
%     centroid_dist_cdf=sort(centroid_dist_mat(:,3)/mean(centroid_dist_mat(:,3)));
%     points_dist_cdf=sort(points_dist_mat(:,3)/mean(points_dist_mat(:,3)));
%     figure(1)
%     hold on;
%     histogram(centroid_dist_mat(:,3)/mean(centroid_dist_mat(:,3)));
%     histogram(points_dist_mat(:,3)/mean(points_dist_mat(:,3))); 
%     figure(2)
%     cdf_centroids_tot=[0,0;centroid_dist_cdf,(1:numel(centroid_dist_cdf))'];
%     semilogy(cdf_centroids_tot(:,1),1-cdf_centroids_tot(:,2)/numel(cdf_centroids_tot(:,2)),'-o','LineWidth',2);
%     hold on
%     cdf_points_tot=[0,0;points_dist_cdf,(1:numel(points_dist_cdf))'];
%     semilogy(cdf_points_tot(:,1),1-cdf_points_tot(:,2)/numel(cdf_centroids_tot(:,2)),'-+','LineWidth',2);
% end
%    