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
% INPUTS: 
%    v_pos=this is an nx2 matrix of vertex positions for the vornoi cell
%    around a point. the v_pos(1,:) is a point at [inf,inf] so we need to
%    remove these points and figure out which vertex needs to replace these
%    infinity points
%
%    connections: This is an Nx1 cell array wherethe ith cell has a list of
%    vertices that sorrounds the ith point in the my_try matrix.
%
%    vx: this is a 2xN matrix of the x locations vertex pairs. Any two
%    points that shart a row are connected.
%
%    vy: this is a 2xN matrix of the y locations vertex pairs. Any two
%    points that shart a row are connected.
%
%    my_try: This is a list of point that we have created the voronoi
%    construction around.
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
%    written by: A. Chieco, UPenn, August  2018
%
%-
function [verts_out,connect_out]=voronoi_vertex_connect(v_pos,connections,vx,vy,p_for_voro)

connect_out=cell(numel(p_for_voro(:,1)),1);
original_id_max=numel(v_pos(:,1));
vert_id_add=original_id_max+1;

v_outside=[0,0,-1];

for i1=1:numel(p_for_voro(:,1))
    old_cons=connections{i1};
    if numel(old_cons(old_cons==1))==0
        connect_out{i1}=old_cons;
    else
        v_incomplete=v_pos(old_cons,:);
        %There should only be one point where this is true
        find_inf=find(old_cons==1);
        if numel(find_inf)==1
            if find_inf(1)==1
                mark1=numel(old_cons);
                mark2=2;               
                new_cons=[0,old_cons(2:end),0];
                v_place1=numel(new_cons);
                v_place2=1;
            elseif find_inf(1)==numel(old_cons)
                mark1=numel(old_cons)-1;                
                mark2=1;
                new_cons=[0,old_cons(1:end-1),0];
                v_place1=numel(new_cons);
                v_place2=1;                
            else
                mark1=find_inf(1)-1;
                mark2=find_inf(1)+1;
                new_cons=[old_cons(1:mark1),0,0,old_cons(mark2:end)];
                v_place1=mark1+1;
                v_place2=mark1+2;
            end
            %we are trying to form a closed loop of vertices from v1 to v2
            for v_find=1:2
                if v_find==1
                    v1=v_pos(old_cons(mark1),:);
                else
                    v1=v_pos(old_cons(mark2),:);
                end
                %we identify which vertices are attached the the current
                %vertex we are investigating. We need to find which of the
                %neighbors lies off of the grid
                d1=sqrt((v1(1)-vx(1,:)).^2+(v1(2)-vy(1,:)).^2);
                d2=sqrt((v1(1)-vx(2,:)).^2+(v1(2)-vy(2,:)).^2);
                v_tag1=find(d1<=1E-12);             
                v_tag2=find(d2<=1E-12);
                if and(numel(v_tag1)==0,numel(v_tag2)==0)
                    v_tag1=find(d1<=1E-8); 
                    v_tag2=find(d2<=1E-8);
                end
                x_same=[vx(1,v_tag1),vx(2,v_tag2)]; y_same=[vy(1,v_tag1),vy(2,v_tag2)];
                if numel(x_same)==0
                    continue
                end
                %we find all the potential new neighbors.
                x_neighbors=[vx(2,v_tag1),vx(1,v_tag2)];
                y_neighbors=[vy(2,v_tag1),vy(1,v_tag2)];
                %the voronoi cell is an equidistant line the seperates two
                %points. We will find which neighbor is in the voronoi cell
                %we are investigating by seeing which neighbor has two
                %distances that are equal, one is the distance to the point
                %whose cell we are trying to complete and the other is the
                %distance to that points neighbor.
                p_dist_mat=zeros(numel(x_neighbors,3));
                for i2=1:numel(x_neighbors)
                    v_already=sqrt((x_neighbors(i2)-v_pos(:,1)).^2+(y_neighbors(i2)-v_pos(:,2)).^2);
                    %This will find the next closest neighbor that is
                    %already identified
                    v_next=sort(v_already(v_already>min(v_already)));
                    v_next=v_next(1);
                    %This is so we don't try and tag the cell with a vertex
                    %id that is already from the list of bulk vertex IDs
                    if or(min(v_already)<1E-10,min(v_already)/v_next<1E-13)
                        p_dist_mat(i2,3)=1E6;
                        continue
                    end
                    dist_points=double(sqrt((x_neighbors(i2)-p_for_voro(:,1)).^2+(y_neighbors(i2)-p_for_voro(:,2)).^2));
                    p_dist_mat(i2,1)=dist_points(i1);
                    dist_points(i1)=2*max(dist_points);
                    %ideally this will be zero for the correct vertex but
                    %for the incorrect neighbors it will not be. Either way
                    %we want to find which other center has the distance
                    %closest to the center we are interrogating.
                    if numel(dist_points(abs(dist_points-p_dist_mat(i2,1))==min(abs(dist_points-p_dist_mat(i2,1)))))>1
                        keyboard                        
                    end
                    p_dist_mat(i2,2)=dist_points(abs(dist_points-p_dist_mat(i2,1))==min(abs(dist_points-p_dist_mat(i2,1))));
                    p_dist_mat(i2,3)=abs(p_dist_mat(i2,1)-p_dist_mat(i2,2));
                end
                %Now we have idenitifed the location of the neighboring
                %vertex. We need to add it to the list of existing vertices
                %if it is not already one there
                neighbor=[x_neighbors(p_dist_mat(:,3)==min(p_dist_mat(:,3))),...
                    y_neighbors(p_dist_mat(:,3)==min(p_dist_mat(:,3)))];
                if numel(neighbor)>2
                        keyboard
                end
                %here we reinitialize v_outside. Now it will only be a
                %matrix of vertex locations and ids that lie outside of the
                %region of interest.
                if v_outside(1,3)==-1
                    v_outside=[neighbor,vert_id_add];
                    v_id_list=vert_id_add;
                    vert_id_add=vert_id_add+1;
                else
                    d_vert_check=sqrt((v_outside(:,1)-neighbor(1)).^2+(v_outside(:,2)-neighbor(2)).^2); 
                    %this means we have found a new vertex. We add it to
                    %the v_outside matrix and add a new ID as well
                    if min(d_vert_check)>1E-12
                        v_outside=[v_outside;neighbor,vert_id_add];
                        v_id_list=vert_id_add;
                        vert_id_add=vert_id_add+1;
                    else
                        v_id_list=v_outside(d_vert_check<1E-12,3);
                    end
                end
%                 Now we have to reconstruct the cell vertices. 
                if v_find==1  
                    new_cons(v_place1)=v_id_list;
                else
                    new_cons(v_place2)=v_id_list;
                end
            end
        else
            keyboard
        end
        if i1==373134
                    keyboard
        end
        connect_out{i1}=new_cons(new_cons>0);
    end
end

verts_out=[[v_pos,(1:original_id_max)'];v_outside];

end

    
