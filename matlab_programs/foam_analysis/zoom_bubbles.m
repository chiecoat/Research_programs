% +
% NAME: watershed_bubbles
%
% PURPOSE:%     
%
% CATEGORY:
%     Foam coarsening
%
% CALLING SEQUENCE:
%    vertex_mat_tot=zoom_vertex_connect(vert_mat)
%
% INPUTS: 
%   verts_info: has columns
%   [id,x,y,theta,vert_type,energy,N_id_1,N_id_2,N_id_3] where N_id
%   is the neighbor ID of -1 if the vertex does not have 3 neighbors.
%
% OPTIONAL INPUTS: (none)
%
% KEYWORD PARAMETERS: (none)
%
% OUTPUTS: 
%       bubbles: a matrix of full bubbles in our foam the columns are 
%
%       verts_out: has columns 
%       [id,x,y,theta,energy,bub_id_1,bub_id_2,bub_id_3,N_id_1,N_id_2,N_id_3]
%       so we remove the "type" of vertex and add the bubble ids a vertex
%       belongs to. If a vertex doesn't belond to 3 bubbles a -1 is used in
%       place of the bubble_id.
% 
% SIDE EFFECTS: (none)
%
% MODIFICATION HISTORY:
%    written by: A. Chieco UPenn, October 2018
%-
function [bubbles,verts_out]=zoom_bubbles(verts_info)

%this initiall seed the bubble matrix
bub_mat=zeros(1,17)-1;
n_verts=numel(verts_info(:,1));
bub_id=1;

for i1=1:n_verts
    vert_id=verts_info(i1,1);
    vert_n_ids_start=verts_info(i1,end-2:end);
    vert_n_ids_start=vert_n_ids_start(vert_n_ids_start>0);
    %I don't think this should ever be zero but this is here just in case.
    %Also is a vertex only hsa one neighbor that means it can't belong to
    %any complete bubble.
    if numel(vert_n_ids_start)<=1
        continue
        %if the vertex has two neighbors then it belongs to one complete bubble
    elseif numel(vert_n_ids_start)==2
        n_bubs=1;
    else
        n_bubs=3;
    end
    %otherwise we try to find a bubble by going around in a loop of
    %vertices. The loop is defined by the direction of the cross product
    %vector and as long as that direction is maintined we remain on one
    %bubble
    for i2=1:n_bubs
        vert_n_ids=vert_n_ids_start;
        %regardless of the vertex neighbors they are all connected through
        %the central vertex we are interrogating.
        if i2==1
            bub_verts=[vert_n_ids(1),vert_id,vert_n_ids(2)];
        elseif i2==2
            bub_verts=[vert_n_ids(1),vert_id,vert_n_ids(3)];
        else
            bub_verts=[vert_n_ids(2),vert_id,vert_n_ids(3)];
        end
        %there are three vertices alrady associated with the bubble
        bub_vert_num=0;
        bub_vert_add=3;
        %we initialize the plane we are on with the first three connected
        %vertices
        vec1=[verts_info(bub_verts(1),2)-verts_info(bub_verts(2),2),...
            verts_info(bub_verts(1),3)-verts_info(bub_verts(2),3)];
        vec2=[verts_info(bub_verts(3),2)-verts_info(bub_verts(2),2),...
            verts_info(bub_verts(3),3)-verts_info(bub_verts(2),3)];
        vec_z=(vec1(1)*vec2(2)-vec1(2)*vec2(1));
        val_z=vec_z/abs(vec_z);
        while bub_vert_num<1
            if bub_vert_add>13
                break
            end
            %we reset our "central" vertex be moving to the end of the
            %bub_verts id list. The id at bub_verts(end-1) is from the
            %direction we came from and and we add ids to check to the end
            %of bub_verts
            %now we reinitialize the the next vertex id and its neighbors.
            %Of the possible neighbors one will be where we came from so we
            %also exlcude that.
            vert_n_ids=verts_info(bub_verts(end),end-2:end);
            vert_n_ids=vert_n_ids(and(vert_n_ids>0,vert_n_ids~=bub_verts(end-1)));
            if numel(vert_n_ids)==0
                bub_vert_num=2;
                break
            end
            %if one other potential neighbor is the starting vertex ID for
            %the bubble then we have completed a loop. That means that we
            %found a complete bubble and can break the while loop. A
            %complete bubble is our "condition 1" and we will add this
            %bubble to a list that needs to be tracked.
            if numel(find(vert_n_ids==bub_verts(1)))==1
                bub_vert_num=1;
            else
                for i3=1:numel(vert_n_ids)
                    bub_verts_check=[bub_verts(end-1:end),vert_n_ids(i3)];
                    vec1=[verts_info(bub_verts_check(1),2)-verts_info(bub_verts_check(2),2),...
                        verts_info(bub_verts_check(1),3)-verts_info(bub_verts_check(2),3)];
                    vec2=[verts_info(bub_verts_check(3),2)-verts_info(bub_verts_check(2),2),...
                        verts_info(bub_verts_check(3),3)-verts_info(bub_verts_check(2),3)];
                    vec_z=(vec1(1)*vec2(2)-vec1(2)*vec2(1));
                    val_check=vec_z/abs(vec_z);
                    if val_check==val_z
                        vert_add=vert_n_ids(i3);
                        bub_vert_add=bub_vert_add+1;
                        break
                    else
                        vert_add=-1;
                    end
                end
                if vert_add>0
                    bub_verts=[bub_verts,vert_add];
                else
                    %Condition 2 means we were unable to close a loop and
                    %make a complete bubble. We want to break the while
                    %loop but also want to indiacte that we do not have a
                    %complete bubble and should not count it as such.
                    bub_vert_num=2;
                end
            end
        end
        if bub_vert_num==1
            bub_vec=zeros(1,17)-1;
            n_sides=numel(bub_verts);
            %before we do any calculations we want to make sure we are not
            %looking at a bubble we have already identified. And we can
            %only possibly have found the same bubble if it has the same
            %number of sides
            bub_mat_check=bub_mat(bub_mat(:,2)==n_sides,:);
            n_check=numel(bub_mat_check);
            if n_check>0
                bub_ids_check=bub_mat_check(:,5:end);
                for i4=1:n_check
                    v_id_match=ismember(bub_verts,bub_ids_check(bub_ids_check>0));
                    %if there is a row where we already have all the same
                    %vertices then that bubble is already found
                    if sum(v_id_match)==n_sides
                        found=1;
                        break
                    else
                        %if not then we have just stumbled upon a bubble
                        %with the same number of sides but at a different
                        %location
                        found=0;
                    end
                end
            else
                found=0;
            end
            if found==1
                continue
            else
                %this is an estimate of the centroid of the bubble based on the area of
                %the polygon made by idneitifiable vertices. This will change but for
                %now is used as a check.
                vert_for_cen=verts_info(bub_verts,2:3);
                n_v_ids=numel(bub_verts);
                a_vec=zeros(1,n_v_ids);
                cx_vec=a_vec;
                cy_vec=a_vec;
                for i5=1:n_v_ids
                    x_i=vert_for_cen(i5,1); y_i=vert_for_cen(i5,2);
                    if i5<n_v_ids
                        x_i_plus1=vert_for_cen(i5+1,1); y_i_plus1=vert_for_cen(i5+1,2);
                    else
                        x_i_plus1=vert_for_cen(1,1); y_i_plus1=vert_for_cen(1,2);
                    end
                    a_vec(i5)=x_i*y_i_plus1-x_i_plus1*y_i;
                    cx_vec(i5)=(x_i+x_i_plus1)*(x_i*y_i_plus1-x_i_plus1*y_i);
                    cy_vec(i5)=(y_i+y_i_plus1)*(x_i*y_i_plus1-x_i_plus1*y_i);
                end
                a_poly=abs(sum(a_vec))/2;
                x_cen=abs(sum(cx_vec))/(6*a_poly);
                y_cen=abs(sum(cy_vec))/(6*a_poly);
                bub_vec(1:5+n_sides-1)=[bub_id,n_sides,x_cen,y_cen,bub_verts];
                bub_id=bub_id+1;
                bub_mat=[bub_mat;bub_vec];
            end
        end
    end
end
%%bubbles outputs a Nx17 row matrix where the columns are 
%[bubble_id, n_sides,x_cen,y_cen,vertex_id1, vertex_id2 ... vertex_id13]
%and a -1 is used as a placeholder for vertex_ids if a bubbles does not
%have 13 sides, which is almost always the case.
bubbles=bub_mat(bub_mat(:,1)>0,:);

%we add three rows  to fill them with the ids of the bubbles that a vertex 
%belongs to
verts_out=[verts_info(:,1:4),verts_info(:,6),zeros(numel(verts_info(:,1)),3)-1,verts_info(:,7:9)];
id_mat=bubbles(:,5:end);
mat_size=size(id_mat);
for i6=1:n_verts
    v_id=i6;
    v_loc=find(id_mat==v_id);
    [bub_row,bub_col]=ind2sub(mat_size,v_loc);   
    if numel(bub_row)>0
        verts_out(i6,6:6+numel(bub_row)-1)=bub_row;
    end    
end

end
