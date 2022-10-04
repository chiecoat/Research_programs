% +
% NAME: square_lattice_make
%
% PURPOSE:
%     This program will solve for h, the "hyperuniformity length" from the variance
%     of a point pattern.
%
% CATEGORY:
%     Hyperuniformity
%
% CALLING SEQUENCE:
%    find_h_length
%
% INPUTS: None at the moment as it is a stand alone program
%
% OPTIONAL INPUTS: (none)
%
% KEYWORD PARAMETERS: (none)
%
% OUTPUTS: Will write a text, that has h computed for different point patterns
%          read in

% SIDE EFFECTS: (none)
%
% MODIFICATION HISTORY:
%    written by: A. Chieco, UPenn, March  2018
%
%-
% clear all
% close all

%Paths where the data will be saved
centroid_path='F:\Chieco\Hyperuniformity\HU foams\voronoi\voro_centroids';
vertex_path='F:\Chieco\Hyperuniformity\HU foams\voronoi\voro_vertices';


num_bubs_vec=[5E3];
rescale=1;
xmin=0; xmax=rescale;
ymin=0; ymax=rescale;

pattern='Einstein_square_delta026';
% 
filename=['F:\Chieco\Hyperuniformity\HU foams\voronoi\points_to_voro\Einstein_square_N5E5 phi_0.79_loop_1 delta_point_0.26 f=0 ad_mag=0.txt'];

for n_count=1:numel(num_bubs_vec)
    for ann=0:0
        cell_txt=num2str(num_bubs_vec(n_count));
        for i0=1:1
%             if ann==0
%                 my_try=rescale*rand(num_bubs_vec(n_count),2);
%             else
%                 my_try=centroid_cells_out(:,1:2);
%             end
            %         a=sobolset(2,'Skip',3E6);
            %         my_try=net(a,num_bubs_vec(n_count));
            %         a=haltonset(2,'Skip',3E6);
            %         my_try=net(a,num_bubs_vec(n_count));
            my_try=dlmread(filename);
            my_try=my_try(:,1:2);
            num_bubs=numel(my_try(:,1));
            %we need to tag the vertice with each other. Thevx and vy will for
            %closed loops but do not give us information about which points lie
            %within the voronoi cell. v_pos and connections does tell us which
            %point is surrounded by which vertices but only has one vertex at
            %infinity,infiinity to describe all of the cell vertices outside
            %some bounded region. We identify the correct loops of vertices
            %around a point.
            [v_pos,connections]=voronoin(my_try);
            [vx,vy]=voronoi(my_try(:,1),my_try(:,2));
            [new_v_pos,new_cells]=voronoi_vertex_connect(v_pos,connections,vx,vy,my_try);
            vert_id_count=numel(new_v_pos(:,1));
            
            
%             figure('pos',[20 50 920 920])
%             ax1 = axes('Position',[0.025 0.025 0.95 0.95]);
%             hold on
%             plot(my_try(:,1),my_try(:,2),'or');
%             plot(v_pos(:,1),v_pos(:,2),'ok')
%             %         plot(new_v_pos(:,1),new_v_pos(:,2),'.m')
%             axis equal
%             xlim([0 rescale])
%             ylim([0 rescale])
%             plot(vx,vy,'b')
%             keyboard
%             close(figure(1))
            %we are going to construct a matrix with cell ID, poisson point and
            %list of ids of the vertices that surround the poisson point. Later
            %we will use these points and vertices to find centroids and areas
            %for the cells
            cells=zeros(num_bubs,20)-1;
            for i1=1:num_bubs
                cells(i1,1:3)=[i1,my_try(i1,:)];
                vert_ids=new_cells{i1};
                v_id_new=-1;
                %Here we are going to adjust the ids of the vertices. We
                %may have to add new vertices to th cell of a bubble
                %because the current ones go over the edge. If that is the
                %cas we will add the new vertex to the vertex list, add an
                %id to the wertex list and switch that new ID into the ids
                %of the vertices for the cell
                for i2=1:numel(vert_ids)
                    if i2<numel(vert_ids)
                        v1=new_v_pos(vert_ids(i2),:);
                        v1_place=i2;
                        v2=new_v_pos(vert_ids(i2+1),:);
                        v2_place=i2+1;
                    else
                        v1=new_v_pos(vert_ids(i2),:);
                        v1_place=i2;
                        v2=new_v_pos(vert_ids(1),:);
                        v2_place=1;
                    end
                    %We need this for a later conditional statement.
                    slope=(v2(2)-v1(2))/(v2(1)-v1(1));
                    b1=v2(2)-slope*v2(1);
                    %this first statement is if at least one of the vertices
                    %lies in the bulk. However there can be points near the
                    %corner that have no vertices in the bulk.
                    if or(and(and(v1(1)>=xmin,v1(1)<=xmax),and(v1(2)>=xmin,v1(2)<=xmax)),...
                            and(and(v2(1)>=xmin,v2(1)<=xmax),and(v2(2)>=xmin,v2(2)<=xmax)))==1
                        vert_find=1;
                        %This means we have a point near the corner that only
                        %has vertices that lie outside the bounding box
                    elseif and(or(v1(1)/v2(1)<0,(v1(1)-1)/(v2(1)-1)<0),or(v1(2)/v2(2)<0,(v1(2)-1)/(v2(2)-1)<0))==1
                        %This is a line that crosses the bottom left corner. Which
                        %means at x=0, y>0 and y<rescale and at y=0 x>0 and x<rescale
                        if and(and(-b1/slope>=0,-b1/slope<=rescale),and(b1>=0,b1<=1))==1
                            vert_find=1;
                            %This is a line that crosses the top left corner. Which
                            %means at x=0, y>0 and y<rescale and at y=rescale x>0 and x<rescale
                        elseif and(and((rescale-b1)/slope>=0,(rescale-b1)/slope<=rescale),and(b1>=0,b1<=1))==1
                            vert_find=1;
                            %This is a line that crosses the top right corner. Which
                            %means at x=rescale, y>0 and y<rescale and at y=rescale x>0 and x<rescale
                        elseif and(and((rescale-b1)/slope>=0,(rescale-b1)/slope<=rescale),and(rescale*slope+b1>=0,rescale*slope+b1<=1))==1
                            vert_find=1;
                            %This is a line that crosses the bottom right corner. Which
                            %means at x=rescale, y>0 and y<rescale and at y=0 x>0 and x<rescale
                        elseif and(and(-b1/slope>=0,-b1/slope<=rescale),and(rescale*slope+b1>=0,rescale*slope+b1<=1))==1
                            vert_find=1;
                            %This means there are two vertices that are connected
                            %across the borders boundaries but not through the
                            %bounding box
                        else
                            vert_find=0;
                        end
                    else
                        vert_find=0;
                    end
                    if vert_find>0
                        %we check the x values of the first vertex
                        if or(v1(1)<0,v1(1)>rescale)==1
                            if v1(1)<0
                                x_new=0;
                            else
                                x_new=rescale;
                            end
                            y_new=slope*x_new+b1;
                            %These check if we are off a diagnoal from the bounding box
                            if y_new>rescale
                                v1(1:2)=[(rescale-b1)/slope,rescale];
                            elseif y_new<0
                                v1(1:2)=[-b1/slope,0];
                            else
                                %This is the standard
                                v1(1:2)=[x_new,y_new];
                            end
                            v_dist=sqrt((v1(1)-new_v_pos(:,1)).^2+(v1(2)-new_v_pos(:,2)).^2);
                            if min(v_dist>1E-12)
                                vert_id_count=vert_id_count+1;
                                v1(3)=vert_id_count;
                                new_v_pos=[new_v_pos;v1];
                            else
                                v1=new_v_pos(v_dist==min(v_dist),:);
                            end
                        end
                        %we check the x values of the second vertex
                        if or(v2(1)<=0,v2(1)>rescale)==1
                            if v2(1)<=0
                                x_new=0;
                            else
                                x_new=rescale;
                            end
                            y_new=slope*x_new+b1;
                            %These check if we are off a diagnoal from the bounding box
                            if y_new>rescale
                                v2(1:2)=[(rescale-b1)/slope,rescale];
                            elseif y_new<0
                                v2(1:2)=[-b1/slope,0];
                            else
                                %This is just the standard
                                v2(1:2)=[x_new,y_new];
                            end
                            v_dist=sqrt((v2(1)-new_v_pos(:,1)).^2+(v2(2)-new_v_pos(:,2)).^2);
                            if min(v_dist>1E-12)
                                vert_id_count=vert_id_count+1;
                                v2(3)=vert_id_count;
                                new_v_pos=[new_v_pos;v2];
                            else
                                v2=new_v_pos(v_dist==min(v_dist),:);
                            end
                        end
                        %We check the y values if they are above the bounding box for
                        %the first vertex
                        if or(v1(2)<=0,v1(2)>rescale)==1
                            if v1(2)<=0
                                y_new=0;
                            else
                                y_new=rescale;
                            end
                            x_new=(y_new-b1)/slope;
                            %We don't need to check off the diagonal because this is
                            %covered in the x check. WE only do the standard value
                            v1(1:2)=[x_new,y_new];
                            v_dist=sqrt((v1(1)-new_v_pos(:,1)).^2+(v1(2)-new_v_pos(:,2)).^2);
                            if min(v_dist>1E-12)
                                vert_id_count=vert_id_count+1;
                                v1(3)=vert_id_count;
                                new_v_pos=[new_v_pos;v1];
                            else
                                v1=new_v_pos(v_dist==min(v_dist),:);
                            end
                        end
                        %We check the y values if they are above the bounding box for
                        %the second vertex
                        if or(v2(2)<=0,v2(2)>rescale)==1
                            if v2(2)<=0
                                y_new=0;
                            else
                                y_new=rescale;
                            end
                            x_new=(y_new-b1)/slope;
                            %We don't need to check off the diagonal because this is
                            %covered in the x check. WE only do the standard value
                            v2(1:2)=[x_new,y_new];
                            v_dist=sqrt((v2(1)-new_v_pos(:,1)).^2+(v2(2)-new_v_pos(:,2)).^2);
                            if min(v_dist>1E-12)
                                vert_id_count=vert_id_count+1;
                                v2(3)=vert_id_count;
                                new_v_pos=[new_v_pos;v2];
                            else
                                v2=new_v_pos(v_dist==min(v_dist),:);
                            end
                        end
                        dist=sqrt((v2(1)-v1(1))^2+(v2(2)-v1(2))^2);
                        if abs(v1(1)-v2(1))>abs(v1(2)-v2(2))
                            xs=(v1(1):(v2(1)-v1(1))/(1E4):v2(1));
                            ys=slope*xs+b1;
                        else
                            ys=(v1(2):(v2(2)-v1(2))/(1E4):v2(2));
                            xs=(ys-b1)/slope;
                        end
                        %                     plot([v1(1);v2(1)],[v1(2),v2(2)],'b','LineWidth',2);
                        %                     plot([v1(1);v2(1)],[v1(2),v2(2)],'og')
                        %                     axis equal
                        %                     xlim([0 rescale])
                        %                     ylim([0 rescale])
                        v_id_new=[v_id_new,v1(3),v2(3)];
                    else
                        v_id_new=[v_id_new,-1,-1];
                    end
                end
                v_id_new=v_id_new(v_id_new>0);
                v_id_new=unique(v_id_new,'stable');
                cells(i1,4:4+numel(v_id_new)-1)=v_id_new;
            end
            %we need to add vertices for the corners and find which are their
            %neighbors.
            id_tot=numel(new_v_pos(:,1));
            corners=[xmin,ymin,id_tot;xmin,ymax,id_tot;xmax,ymax,id_tot;xmax,ymin,id_tot];
            corner_ids=zeros(4,3);
            id_list=cells(:,4:end);
            corner_id_add=1;
            for i3=1:4
                corner_xy=new_v_pos(and(new_v_pos(:,1)==corners(i3,1),...
                    new_v_pos(:,2)==corners(i3,2)),:);
                if numel(corner_xy)>0
                    corners(i3,:)=[-1,-1,-1];
                    continue
                end
                corners(i3,3)=corners(i3,3)+corner_id_add;
                corner_id_add=corner_id_add+1;
                x_con_tot=new_v_pos(new_v_pos(:,2)==corners(i3,2),:);
                x_look=abs(x_con_tot(:,1)-corners(i3,1));
                x_keep=x_con_tot(x_look==min(x_look),:);
                y_con_tot=new_v_pos(new_v_pos(:,1)==corners(i3,1),:);
                y_look=abs(y_con_tot(:,2)-corners(i3,2));
                y_keep=y_con_tot(y_look==min(y_look),:);
                corner_ids(i3,:)=[corners(i3,3),x_keep(3),y_keep(3)];
                id_rows1=find(id_list==x_keep(3));
                [ind_row1,ind_col1]=ind2sub(size(id_list),id_rows1);
                id_rows2=find(id_list==y_keep(3));
                [ind_row2,ind_col2]=ind2sub(size(id_list),id_rows2);
                %We only care about the row since that won't change. ONe of the
                %rows should be the same since that is the cell the two
                %vertices share
                bub_corner=intersect(ind_row1,ind_row2);
                cell_tot=cells(bub_corner,:);
                if numel(cell_tot(:,1))>1
                    for i4=1:numel(cell_tot(:,1))
                        cell_ids_in=cell_tot(i4,4:end);
                        %This means we have two vertices that intersct the
                        %boundary. The cell that only has two vertices in not a
                        %shape so it needs the corner to complete one. This is
                        %where we add the corner point.
                        if numel(new_v_pos(cell_ids_in(cell_ids_in>0),1))==2
                            bub_corner=bub_corner(i4);
                            verts_in_cell=cell_ids_in(1:2);
                            new_verts_in_cell=[verts_in_cell,corners(i3,3)];
                            break
                        else
                            continue
                        end
                    end
                else
                    cell_ids=cell_tot(4:end);
                    verts_in_cell=cell_ids(cell_ids>0);
                    spots=[find(cell_ids==x_keep(3)),find(cell_ids==y_keep(3))];
                    spot1=min(spots);
                    spot2=max(spots);
                    %we need to know where to put the corner vertex.
                    if and(spot1==1,spot2==numel(verts_in_cell))==1
                        new_verts_in_cell=[verts_in_cell,corners(i3,3)];
                    elseif and(spot1==1,spot2==2)==1
                        new_verts_in_cell=[verts_in_cell(1),corners(i3,3),verts_in_cell(2:end)];
                    else
                        new_verts_in_cell=[verts_in_cell(1:spot1),corners(i3,3),verts_in_cell(spot2:end)];
                    end
                end
                cells(bub_corner,4:4+numel(new_verts_in_cell)-1)=new_verts_in_cell;
            end
            new_v_pos=[new_v_pos;corners(corners(:,3)>0,:)];
            %Now we need to tag the vertices and only keep the ones that lie
            %within our bounding box.
            new_v_pos=[new_v_pos,zeros(numel(new_v_pos(:,1)),1)-1];
            new_v_pos(and(and(new_v_pos(:,1)>=0,new_v_pos(:,1)<=rescale),and(new_v_pos(:,2)>=0,new_v_pos(:,2)<=rescale)),4)=...
                (1:numel(new_v_pos(and(and(new_v_pos(:,1)>=0,new_v_pos(:,1)<=rescale),and(new_v_pos(:,2)>=0,new_v_pos(:,2)<=rescale)),4)))';
            %Now we loop through each cell (1) retag the vertices and (2) find
            %the area of the cell. We will also find the centroid of the
            %voronoi cell and make a seperate list of centroids with cell area
            %these matrices will just hold x,y,area,n_sides for patterns made by either
            %the poisson pattern or by the centroid pattern
            points_cells_out=zeros(numel(cells(:,1)),4);
            centroid_cells_out=zeros(numel(cells(:,1)),4);
            %we are adding a row for the cell area, the number of sides, and the
            %centroid of the cell (as opposed to the poisson point the cell is
            %generated around)
            cells_tot_out=[cells,zeros(numel(cells(:,1)),4)-1];
            for i5=1:numel(cells(:,1))
                cell_ids_old=cells(i5,4:end);
                verts_def=new_v_pos(cell_ids_old(cell_ids_old>0),:);
                cell_ids_new=verts_def(:,4);
                n_sides=numel(cell_ids_new);
                %Now we find the area and centroid of the cell
                area_vec=zeros(1,n_sides);
                cen_x_vec=zeros(1,n_sides);
                cen_y_vec=zeros(1,n_sides);
                for i6=1:n_sides
                    if i6<n_sides
                        v_xy1=verts_def(i6,1:2);
                        v_xy2=verts_def(i6+1,1:2);
                    else
                        v_xy1=verts_def(i6,1:2);
                        v_xy2=verts_def(1,1:2);
                    end
                    cross_terms=v_xy1(1)*v_xy2(2)-v_xy2(1)*v_xy1(2);
                    area_vec(i6)=cross_terms;
                    cen_x_vec(i6)=(v_xy1(1)+v_xy2(1))*cross_terms;
                    cen_y_vec(i6)=(v_xy1(2)+v_xy2(2))*cross_terms;
                end
                %The area of the polygon is signed depending on if the vertices
                %go clockwise or counterclockwise. We do not care about the
                %order so we take the absolute value
                area_tot=sum(area_vec)/2;
                cen_x=sum(cen_x_vec)/(6*area_tot);
                cen_y=sum(cen_y_vec)/(6*area_tot);
                points_cells_out(i5,:)=[cells(i5,2:3),abs(area_tot),n_sides];
                centroid_cells_out(i5,:)=[cen_x,cen_y,abs(area_tot),n_sides];
                cells_tot_out(i5,1:7)=[i5,n_sides,abs(area_tot),cells(i5,2:3),cen_x,cen_y];
                cells_tot_out(i5,8:8+n_sides-1)=cell_ids_new;
            end
            if abs(sum(centroid_cells_out(:,3))-1)>1E-8
                keyboard
            else
                abs(sum(centroid_cells_out(:,3))-1)
                %there are four files we want to write out, one for the whole
                %list of information for the cells, one has the ids,x,y
                %positions of the vertices, one is the poisson pattern and the
                %corresponding cell areas and the last one is the centroid
                %patter and the corresponding area
                centroid_info=[centroid_path '\x_y_area_N lists\' pattern ' centroid_xyAN N_cells_' cell_txt ' run_' num2str(i0) ' .txt'];
                dlmwrite(centroid_info,centroid_cells_out,'Delimiter',',','Precision','%1.15e','newline','pc');
                points_info=[centroid_path '\x_y_area_N lists\' pattern ' pointsINI_xyAN N_cells_' cell_txt ' run_' num2str(i0) ' .txt'];
                dlmwrite(points_info,points_cells_out,'Delimiter',',','Precision','%1.15e','newline','pc');
                cell_info_tot=[centroid_path '\all info list\' pattern ' ID_NA_xyp_xyc_vertID N_cells_' cell_txt ' run_' num2str(i0) ' .txt'];
                dlmwrite(cell_info_tot,cells_tot_out,'Delimiter',',','Precision','%1.15e','newline','pc');
                verts_out=[new_v_pos(new_v_pos(:,3)>0,1:2),new_v_pos(new_v_pos(:,3)>0,4)];
                vert_info_tot=[vertex_path '\' pattern ' vertex_xyID N_cells_' cell_txt ' run_' num2str(i0) '.txt'];
                dlmwrite(vert_info_tot,verts_out,'Delimiter',',','Precision','%1.15e','newline','pc');
            end
        end
    end
end


