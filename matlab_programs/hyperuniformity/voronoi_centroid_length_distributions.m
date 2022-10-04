readPath='F:\Chieco\Hyperuniformity\HU foams\voronoi\';
pattern={'Halton','Einstein_square_delta025','Einstein_square','poisson'};


for i1=1:numel(pattern)
    centroid_dist_mat=[-1,-1,-1];
    points_dist_mat=[-1,-1,-1];
    pat=pattern{i1};
    cell_mat=dlmread([readPath '\voro_centroids\all info list\' pat ' ID_NA_xyp_xyc_vertID N_cells_5000 run_1 .txt']);
    n_bubs=numel(cell_mat(:,1));
    vert_mat=cell_mat(:,8:end);
    v_sz=size(vert_mat);
    for i2=1:n_bubs
        id=cell_mat(i2,1);
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
            if numel(centroid_dist_mat(and(centroid_dist_mat(:,1)==id_pair(1),centroid_dist_mat(:,2)==id_pair(2)),:))>1
               continue
            end
            dist_cens=sqrt((cell_mat(id_pair(1),4)-cell_mat(id_pair(2),4))^2+...
                           (cell_mat(id_pair(1),5)-cell_mat(id_pair(2),5))^2);
            dist_points=sqrt((cell_mat(id_pair(1),6)-cell_mat(id_pair(2),6))^2+...
                           (cell_mat(id_pair(1),7)-cell_mat(id_pair(2),7))^2);
            centroid_dist_mat=[centroid_dist_mat;id_pair,dist_cens];
            points_dist_mat=[points_dist_mat;id_pair,dist_points];
        end
    end
    centroid_dist_mat=centroid_dist_mat(2:end,:);
    points_dist_mat=points_dist_mat(2:end,:);
    figure(5)
    hold on;
    histogram(centroid_dist_mat(:,3)/mean(centroid_dist_mat(:,3)));
    histogram(points_dist_mat(:,3)/mean(points_dist_mat(:,3))); 
end
    
    
    


