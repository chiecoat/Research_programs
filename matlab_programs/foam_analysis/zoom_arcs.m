% +
% NAME: watershed_bubbles
%
% PURPOSE:%     
%
% CATEGORY:
%     Foam coarsening
%
% CALLING SEQUENCE:
%    bubbles=watershed_bubble(vertex_info)
%
% INPUTS: 
%    verts_info: this is a Nx8 matrix where the columns are 
%    [id_vert,x_vert,y_vert,theta,energy,cell_id_1,cell_id_2,cell_id_3] 
%
% OPTIONAL INPUTS: (none)
%
% KEYWORD PARAMETERS: (none)
%
% OUTPUTS: Will write files to destination folders
% 
% SIDE EFFECTS: (none)
%
% MODIFICATION HISTORY:
%    written by: A. Chieco UPenn, October 2018
%-
function [arcs,arcs_err,arc_imi]=zoom_arcs(foam_imi,vert_mat,vert_errors,dist_thresh)

%each bubble is already sorted in a way so that the vertices are
%counter-clockwise so pairs of adjacent vertex IDs are actual vertex
%neighbors
arc_imi=foam_imi;
arc_imi_out=arc_imi*0;
L_sys_x=numel(foam_imi(1,:));
L_sys_y=numel(foam_imi(:,1));
arcs=zeros(1,5);
arcs_err=zeros(1,5);

% imshow(foam_imi/255)
% hold on

for i1=1:numel(vert_mat(:,1))
      v_id_cen=vert_mat(i1,1);
      v_id_neighbors=vert_mat(i1,end-2:end);
      vert_ids=v_id_neighbors(v_id_neighbors>0);
    for i2=1:numel(vert_ids)
        %We identify the two vertices we are trying to connect
         p1_xy=vert_mat(v_id_cen,2:3);
         p1_xy_err=vert_errors(v_id_cen,1:2);
         p2_xy=vert_mat(vert_ids(i2),2:3);
         p2_xy_err=vert_errors(vert_ids(i2),1:2);
         vec_ids=sort([v_id_cen,vert_ids(i2)]);
        %Each vertex pair appears twice as we perform this code, because
        %pairs of vertices make the face of two bubbles.
        %Below we make sure we do not double count the vertex pairs
        if i1>1
            my_check=arcs(and(arcs(:,1)==vec_ids(:,1),...
                arcs(:,2)==vec_ids(:,2)),:);
            if numel(my_check)>0
                continue
            end
        end   
        %These are the arc equations we want to fit the actual films. We
        %need to know the two points that are the vertex locations 
        [arc_eqn,x_vals,y_vals]=zoom_arc_minimum_find(p1_xy,p2_xy,foam_imi,100);
        if and(sum(isnan(p1_xy_err))>0,sum(isnan(p2_xy_err))==0)==1
            p1_xy_err=p2_xy_err;
        elseif and(sum(isnan(p1_xy_err))==0,sum(isnan(p2_xy_err))>0)==1
            p2_xy_err=p1_xy_err;
        end
        %Now we also fit for the arcs that make the error        
        [arc_eqn_err_mat,xy_pinch,xy_pull]=zoom_arc_minimum_error_find(p1_xy,p2_xy,p1_xy_err,p2_xy_err,foam_imi,100);
        arc_eqn_err=arc_eqn_err_mat(1,:);
        for i3=1:2
            if i3==2
                xy_vals=[x_vals,y_vals];
                conv_mat_val=1;
            else
                xy_vals=xy_pinch;
                conv_mat_val=2;
            end
            circ_xy_tot=[round(xy_vals(:,1)),round(xy_vals(:,2))];
            circ_xs_keep=circ_xy_tot(and(circ_xy_tot(:,1)>0,circ_xy_tot(:,1)<L_sys_x),:);
            if numel(circ_xs_keep)<2
                continue
            end
            circ_xy_ind=circ_xs_keep(and(circ_xs_keep(:,2)>0,circ_xs_keep(:,2)<L_sys_y),:);
            arc_ind=sub2ind([L_sys_y,L_sys_x],circ_xy_ind(:,2),circ_xy_ind(:,1));
            arc_imi_out(arc_ind)=1;
            insert_mat=conv2(arc_imi_out,zeros(9,9)+conv_mat_val,'same');
            arc_imi(insert_mat>0)=conv_mat_val;
            arc_imi_out=arc_imi_out*0;  
%             plot(circ_xy_ind(:,1),circ_xy_ind(:,2),'or')
        end
        %our arc matrix will have x_0,y_0,R for the circlew we fit to the
        %films follwoing the vertex ids
        arcs=[arcs;[vec_ids,arc_eqn]];
        %The arc_err matrix will have [x_pinch,y_pinch,R_pinch, x_pull,y_pull,R_pull]
        %following the vertex ids
        arcs_err=[arcs_err;[vec_ids,arc_eqn_err]];
    end
end

arcs=arcs(2:end,:);
arcs_err=arcs_err(2:end,:);


end


