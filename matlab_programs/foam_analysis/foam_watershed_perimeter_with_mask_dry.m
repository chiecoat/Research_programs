clear all
close all


%This is the folder we read the dry data from
date_vec={'2022-09-26 nitrogen too mono'};
loop_start_vec=[1944];
loop_end_vec=[2250];
% loop_start_vec=[99,99,99];
% loop_end_vec=[100,100,100];

side_color_mat_1=jet(21);
% side_color_mat_1=circshift(side_color_mat_1,[3 0]);
% insert=side_color_mat_1(4,:);
% side_color_mat_1(4,:)=side_color_mat_1(end-1,:);
% side_color_mat_1(end-1,:)=insert;

% side_color_mat_2=round(jet(20));
% side_color_mat_2=circshift(side_color_mat_2,[8 0]);

side_color_mat=[side_color_mat_1((1:5:21),:);side_color_mat_1((2:5:17),:);...
                side_color_mat_1((3:5:18),:);side_color_mat_1((4:5:19),:);...
                side_color_mat_1((5:5:20),:)];


for dates=1:1 
    folder_titles=date_vec{dates};
    %First we read in all of the files with the bubbles information
    imiPath=['E:\foam coarsening\non-scaling state\' folder_titles '\binary full'];
    filePath=['E:\Chieco\non-scaling state\' folder_titles '\bubble data'];
    
    num_start=0;
    loop_start=loop_start_vec(dates);
    loop_end=loop_end_vec(dates);
    % loop_end=1089;
    
    for num_foam=loop_start:6:loop_end
        foam_mask=imread(['E:\foam coarsening\non-scaling state\' folder_titles '\mask.png']);
        foam_mask_points=find(foam_mask==0);
        clear foam_mask
        foam_mask_edge=imread(['E:\foam coarsening\non-scaling state\' folder_titles '\mask_edge.png']);

        extra=20;
        foam_mask_edge(extra+1:end-5*extra-1,extra+1)=0;
        foam_mask_edge(extra+1:end-5*extra-1,end-extra-1)=0;
        foam_mask_edge(extra+1,extra+1:end-extra-1)=0;
        foam_mask_edge(end-5*extra-1,extra+1:end-extra-1)=0;
        foam_mask_edge(foam_mask_points)=255;
        foam_mask_edge(1:extra,:)=255;
        foam_mask_edge(end-5*extra:end,:)=255;
        foam_mask_edge(:,1:extra)=255;
        foam_mask_edge(:,end-extra:end)=255;        
        foam_mask_edge_points=find(foam_mask_edge==0);
        clear foam_mask_edge        
        
        if num_foam<10
            vertexPath=[imiPath '\binary DSC_000' num2str(num_start+num_foam) '.png'];
        elseif and(num_foam>=10,num_foam<100)==1
            vertexPath=[imiPath '\binary DSC_00' num2str(num_start+num_foam) '.png'];
        elseif and(num_foam>=100,num_foam<1000)==1
            vertexPath=[imiPath '\binary DSC_0' num2str(num_start+num_foam) '.png'];
        else
            vertexPath=[imiPath '\binary DSC_' num2str(num_start+num_foam) '.png'];
        end
        vertex_imi=imread(vertexPath);
        vertex_imi(foam_mask_points)=0;
        vertex_imi(foam_mask_edge_points)=255;        
        L_sys_x=numel(vertex_imi(1,:));
        L_sys_y=numel(vertex_imi(:,1));
        
        vertex_imi(vertex_imi<150)=0;
        vertex_imi(vertex_imi>150)=255;        
        
        vertex_imi(1:extra,:)=0;
        vertex_imi(end-5*extra:end,:)=0;
        vertex_imi(:,1:extra)=0;
        vertex_imi(:,end-extra:end)=0;
        vertex_imi(extra+1:end-5*extra-1,extra+1)=255;
        vertex_imi(extra+1:end-5*extra-1,end-extra-1)=255;
        vertex_imi(extra+1,extra+1:end-extra-1)=255;
        vertex_imi(end-5*extra-1,extra+1:end-extra-1)=255;
%         %we are going to pad watershed imi and then draw a box around the foam
        watershed_full=vertex_imi;
        
        Ld_full=watershed(watershed_full);
        Ld_full(Ld_full==Ld_full(1,end))=1;
                
        Ld_full(extra+1:end-5*extra-1,extra+1)=0;
        Ld_full(extra+1:end-5*extra-1,end-extra-1)=0;
        Ld_full(extra+1,extra+1:end-extra-1)=0;
        Ld_full(end-5*extra-1,extra+1:end-extra-1)=0;
        
        Ld_full(foam_mask_points)=1;
        Ld_full(foam_mask_edge_points)=0;  
        
        Ld_full(1:extra,:)=1;
        Ld_full(end-5*extra:end,:)=1;
        Ld_full(:,1:extra)=1;
        Ld_full(:,end-extra:end)=1;

%         rgb = label2rgb(Ld_full,'jet','k','shuffle');
%         imshow(rgb);
%         hold on
        
        stats=regionprops(Ld_full,'Area','Centroid','Circularity','Eccentricity','Perimeter');
        %We made a cell that represents the perimeter of the foam. We elimiate
        %that cell which is always the first one.
        areas=[stats.Area];
        n_regions=numel(areas);
        circ=[stats.Circularity];
        ecc=[stats.Eccentricity];
        perim=[stats.Perimeter];
        %we need to identify the particle centers
        cens=[stats.Centroid];
        x_cens=cens((1:2:2*n_regions-1));
        y_cens=cens((2:2:2*n_regions));
        
        %we find the x and y of foam frame, plot the x and y in a figure
        thresh = find(Ld_full==0);
        [y_vec,x_vec]=ind2sub([L_sys_y,L_sys_x],thresh);
        xy_mat=[x_vec,y_vec];
        
        [y_edges,x_edges]=ind2sub([L_sys_y,L_sys_x],foam_mask_edge_points);
        xy_edge_mat=[x_edges,y_edges];
        
        clear foam_mask_edge_points
        clear foam_mask_points
        
        %here we find the location of the foam vertices
        xy_tag=[xy_mat,zeros(numel(xy_mat(:,1)),3)];
        xy_tag(:,5)=-1;
        for i2=1:numel(xy_mat(:,1))
            % establish a small box named Ld_small at each x and y
            if or(xy_mat(i2,2)-1<1,xy_mat(i2,2)+1>numel(Ld_full(:,1)))==1
                continue
            elseif or(xy_mat(i2,1)-2<1,xy_mat(i2,1)+1>numel(Ld_full(1,:)))==1
                continue
            end
            Ld_small=Ld_full((xy_mat(i2,2)-1):(xy_mat(i2,2)+1),(xy_mat(i2,1)-1):(xy_mat(i2,1)+1));
            % find all the different values in the small box. Unique returns the
            % values sorted
            unique_vals=unique(Ld_small,'sorted');
            % find the non-zero value
            bub_id=find(unique_vals~=0);            
            if numel(bub_id)~=3
                continue
            end
            a_check=areas(unique_vals(bub_id));
            if numel(a_check(a_check<10))>0
                continue
            end
            % add the non-zero values into third through n columns of i2 row
            xy_tag(i2,3:3+numel(bub_id)-1)=unique_vals(bub_id);
        end
        % keep the first and second columns of xy_tag if the value of fifth column is greater than
        % zero
        xy_keep=xy_tag(xy_tag(:,5)>0,:);
%         hold on
%         plot(xy_keep(:,1),xy_keep(:,2),'ok')
        
        %here we the bubble ids thatoverlap with the edge of the mask
%         verts_edge=xy_tag(xy_tag(:,3)==1,3:5);
%         bub_edge_ids=reshape(verts_edge,[numel(verts_edge),1]);
%         bubs_elim=unique(bub_edge_ids);

        bubs_elim=[-1];
        for i2=1:numel(xy_edge_mat(:,1))
            % establish a small box named Ld_small at each x and y
            Ld_small=Ld_full((xy_edge_mat(i2,2)-1):(xy_edge_mat(i2,2)+1),(xy_edge_mat(i2,1)-1):(xy_edge_mat(i2,1)+1));
            % find all the different values in the small box. Unique returns the
            % values sorted
            unique_vals=unique(Ld_small,'sorted');
            % find the non-zero value
            bub_id=find(unique_vals~=0);
            % add the non-zero values into third through n columns of i2 row
            bubs_elim=[bubs_elim;unique_vals(bub_id)];
        end
        bubs_elim=bubs_elim(2:end);
        bubs_elim=unique(bubs_elim);
        
        n_bubs=numel(areas);
        %we put all the data tgether in to one matrix
        area_write=[x_cens',y_cens',zeros(n_bubs,1),areas',circ',ecc',perim'];
        
        bub_ids_list=[0,0,0];
        for i3=1:numel(areas)
            if numel(bubs_elim(bubs_elim==i3))==1
                area_write(i3,3)=-1;
            else
                sides1=find(xy_keep(:,3)==i3);
                sides2=find(xy_keep(:,4)==i3);
                sides3=find(xy_keep(:,5)==i3);
                %We know how many sides there are for the bubbles
                area_write(i3,3)=numel(sides1)+numel(sides2)+numel(sides3);
                length_add=ceil(sqrt(area_write(4)/pi)/10);
                perim_check=[sides1;sides2;sides3];
                if numel(perim_check)<3
                    continue
                end
                vert_locs=xy_keep(perim_check,:);
                if or(max(vert_locs(:,1))-min(vert_locs(:,1))==0,max(vert_locs(:,2))-min(vert_locs(:,2))==0)==1
                    continue
                end
                perim_sort=convhull(vert_locs(:,1),vert_locs(:,2));
                perim_values=zeros(1,area_write(i3,3));
                for i4=1:numel(perim_sort)-1
                    vert1=vert_locs(perim_sort(i4),:);
                    vert2=vert_locs(perim_sort(i4+1),:);
                    %this finds the two common vector elements in vert1 and
                    %vert2. The two common values are the ids of the bubbles
                    %that the film connects. The elements are returned sorted
                    bub_ids=intersect(vert1(3:5),vert2(3:5));
                    if numel(bub_ids)~=2
                        continue
                    end
                    if and(bub_ids_list(:,1)==bub_ids(1),bub_ids_list(:,2)==bub_ids(2))==1
                        keyboard
                        perim_values(i4)=bub_ids_list(and(bub_ids_list(:,1)==bub_ids(1),bub_ids_list(:,2)==bub_ids(2)),3);
                    else
                        if abs(vert1(1)-vert2(1))>abs(vert1(2)-vert2(2))
                            film_box=Ld_full(min([vert1(2),vert2(2)])-length_add:max([vert1(2),vert2(2)])+length_add,...
                                min([vert1(1),vert2(1)]):max([vert1(1),vert2(1)]));
                            film_between=zeros(1,numel(film_box(1,:)));
                            for i5=1:numel(film_box(1,:))
                                vals1=find(film_box(:,i5)==bub_ids(1));
                                vals2=find(film_box(:,i5)==bub_ids(2));
                                if or(numel(vals1)==0,numel(vals2)==0)==1
                                    continue
                                else
                                    dist_1=abs(max(vals1)-min(vals2));
                                    dist_2=abs(max(vals2)-min(vals1));
                                end
                                film_between(i5)=min([dist_1,dist_2])-1;
                            end
                            perim_values(i4)=sum(film_between);
                            bub_ids_list=[bub_ids_list;bub_ids,sum(film_between)];
                        else
                            film_box=Ld_full(min([vert1(2),vert2(2)]):max([vert1(2),vert2(2)]),...
                                min([vert1(1),vert2(1)])-length_add:max([vert1(1),vert2(1)])+length_add);
                            film_between=zeros(1,numel(film_box(:,1)));
                            for i5=1:numel(film_box(:,1))
                                vals1=find(film_box(i5,:)==bub_ids(1));
                                vals2=find(film_box(i5,:)==bub_ids(2));
                                if or(numel(vals1)==0,numel(vals2)==0)==1
                                    continue
                                else
                                    dist_1=abs(max(vals1)-min(vals2));
                                    dist_2=abs(max(vals2)-min(vals1));
                                end
                                film_between(i5)=min([dist_1,dist_2])-1;
                            end
                            perim_values(i4)=sum(film_between);
                            bub_ids_list=[bub_ids_list;bub_ids,sum(film_between)];
                        end
                    end
                end
                area_write(i3,4)=area_write(i3,4)+sum(perim_values)/2+area_write(i3,3)/3;
            end
        end
        
        area_keep=area_write(area_write(:,3)>2,:);
        
        bubPath=[filePath '\bubbles imi_' num2str(num_start+num_foam) '_xyNA_CirEccPer_data_full.txt'];
        dlmwrite(bubPath,area_keep,'Delimiter',',','Precision','%1.8e','newline','pc');
        
        n_side_color_coordinate=zeros(numel(Ld_full(:,1)),numel(Ld_full(1,:)),3);
        n_side_label_r=n_side_color_coordinate(:,:,1);
        n_side_label_g=n_side_color_coordinate(:,:,2);
        n_side_label_b=n_side_color_coordinate(:,:,3);
        for label_bub=1:numel(area_keep(:,1))
            label=Ld_full(floor(area_keep(label_bub,2)),floor(area_keep(label_bub,1)));
            n_side_label_r(Ld_full==label)=side_color_mat(area_keep(label_bub,3)-2,1);
            n_side_label_g(Ld_full==label)=side_color_mat(area_keep(label_bub,3)-2,2);
            n_side_label_b(Ld_full==label)=side_color_mat(area_keep(label_bub,3)-2,3);
        end
%         figure(2)        
        n_side_color_coordinate(:,:,1)=n_side_label_r;
        n_side_color_coordinate(:,:,2)=n_side_label_g;
        n_side_color_coordinate(:,:,3)=n_side_label_b;
%         imshow(n_side_color_coordinate);
         
        imwrite(n_side_color_coordinate,side_color_mat_1,['E:\Chieco\non-scaling state\' folder_titles '\watershed images\full watershed imi_' num2str(num_start+num_foam) '.png'])
         
%         pause(45)
    end
end
