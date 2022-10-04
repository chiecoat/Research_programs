% +
% NAME: foam_watershed_2d
%
% PURPOSE:%     
%
% CATEGORY:
%     Foam coarsening
%
% CALLING SEQUENCE:
%    foam_reconstruct_2d
%
% INPUTS: 
%    p1_xy: This is the x and y positions of the vertex we are trying to
%    fit circular arc to.
%
%    p2_xy: This is the other pair if x and y positions of the vertex we 
%    are trying to fit circular arc to.
%
%    error: This is a 2x2 matrix of errirs in the x,y locations of 
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
%    written by: A. Chieco, UPenn, September 2017
%-
function [arc_eq_out,xy_pinch,xy_pull]=zoom_arc_minimum_error_find(p1_xy,p2_xy,p1_err,p2_err,foam_imi,min_thresh)

L_sys_x=numel(foam_imi(1,:));
L_sys_y=numel(foam_imi(:,1));

%First we move the x positions of each vertex around. We either "pinch" or
%"pull" where we move the x positions closer together or further apart
if p1_xy(1)>p2_xy(1)
    p1_x_pinch=p1_xy(1)-p1_err(1);
    p2_x_pinch=p2_xy(1)+p2_err(1);
    p1_x_pull=p1_xy(1)+p1_err(1);
    p2_x_pull=p2_xy(1)-p2_err(1);
else
    p1_x_pinch=p1_xy(1)+p1_err(1);
    p2_x_pinch=p2_xy(1)-p2_err(1);
    p1_x_pull=p1_xy(1)-p1_err(1);
    p2_x_pull=p2_xy(1)+p2_err(1);
end
%We do the same thing for the y positions
if p1_xy(2)>p2_xy(2)
    p1_y_pinch=p1_xy(2)-p1_err(2);
    p2_y_pinch=p2_xy(2)+p2_err(2);
    p1_y_pull=p1_xy(2)+p1_err(2);
    p2_y_pull=p2_xy(2)-p2_err(2);
else
    p1_y_pinch=p1_xy(2)+p1_err(2);
    p2_y_pinch=p2_xy(2)-p2_err(2);
    p1_y_pull=p1_xy(2)-p1_err(2);
    p2_y_pull=p2_xy(2)+p2_err(2);
end

arc_eq_out=zeros(2,3);
for errs=1:2
    if errs==1
        del_x=p1_err(1)-p2_err(1);
        del_y=p1_err(2)-p2_err(2);
    else
        del_x=-(p1_err(1)-p2_err(1));
        del_y=-(p1_err(2)-p2_err(2));
    end
    %we need to know the distance between the two vertices
    d_verts=sqrt((p1_xy(1)-p2_xy(1))^2+(p1_xy(2)-p2_xy(2))^2);
    %We find the midpoint between the two vertices
    pm_xy=[(p1_xy(1)+p2_xy(1))/2,(p1_xy(2)+p2_xy(2))/2];
    slope_m=(p2_xy(2)-p1_xy(2))/(p2_xy(1)-p1_xy(1));
    %Slope 0 will give us an infinite radius, so we consider vertices that lie
    %on in the same plane
    if slope_m==0
        slope_m=1E-12;
    end
    %We find the line that perpendicularly bisects the two vertices
    slope_bis=-1/slope_m;
    y_int_bis=pm_xy(2)-slope_bis*pm_xy(1);
    denom=sqrt(1+slope_bis^2);
    %now we scan the image for points we want to try and fit circles to. The
    %vertices are either "horizontally oriented" or "vertically oriented" and
    %we have to distinguish between the two cases
    xs_test=abs(p1_xy(1)-p2_xy(1));
    ys_test=abs(p1_xy(2)-p2_xy(2));
    %"Horizontally oriented" vertices lie further apart in x than y
    if xs_test>=ys_test
        %These two points make a line that bisects the the vertices
        x_seek_plus=pm_xy(1)+d_verts/(2*denom);
        x_seek_minus=pm_xy(1)-d_verts/(2*denom);
        %There is a chance the two points are colinear
        if x_seek_plus~=x_seek_minus
            xs_seek=(x_seek_plus:(x_seek_minus-x_seek_plus)/d_verts:x_seek_minus);
            ys_seek=slope_bis*xs_seek+y_int_bis;
        else
            ys_seek=(-d_verts/2:1:d_verts/2)+pm_xy(2);
            xs_seek=zeros(1,numel(ys_seek))+x_seek_plus;
        end
    else
        %"Vertically oriented" vertices lie further apart in y than x
        %These two points make a line that bisects the the vertices
        y_seek_plus=pm_xy(2)+slope_bis*d_verts/(2*denom);
        y_seek_minus=pm_xy(2)-slope_bis*d_verts/(2*denom);
        %There is a chance the two points are colinear
        if y_seek_plus~=y_seek_minus
            ys_seek=(y_seek_plus:(y_seek_minus-y_seek_plus)/d_verts:y_seek_minus);
            xs_seek=(ys_seek-y_int_bis)/slope_bis;
        else
            xs_seek=(-d_verts/2:1:d_verts/2)+pm_xy(1);
            ys_seek=zeros(1,numel(xs_seek))+y_seek_plus;
        end
    end
    %We convert indices (x,y) positions [which are also column,row indices]
    %into a single subscript to interrogate the matrix
    mat_index_round=[round(xs_seek'),round(ys_seek')];
    mat_index_xs=mat_index_round(and(mat_index_round(:,1)>0,mat_index_round(:,1)<=L_sys_x),:);
    mat_index_tot=mat_index_xs(and(mat_index_xs(:,2)>0,mat_index_xs(:,2)<=L_sys_y),:);
    val_ind=sub2ind([L_sys_y,L_sys_x],mat_index_tot(:,2),mat_index_tot(:,1));
    my_vals=foam_imi(val_ind);
    ids=(1:numel(my_vals));
    window=min(numel(my_vals)/10,20);
    vals_smooth=smoothdata(my_vals,'gaussian',window);
    %now we identify the region between the films that we want to fit a circle
    %to. To do this we sucessively cut away regions based on known patterns in
    %brightness starting with the first local minima. This allows us to remove
    %the grayscale values that are associated with the interiro of the bubbles
    [minima]=islocalmin(vals_smooth);
    first_pass_min=[ids(minima)',vals_smooth(minima)];
    my_keep_min=first_pass_min(first_pass_min(:,2)<=min_thresh,:);
    
    ids_pass=(my_keep_min(1,1):my_keep_min(end,1));
    vals_smooth_pass=vals_smooth(my_keep_min(1,1):my_keep_min(end,1));
    
    [maxima]=islocalmax(vals_smooth_pass);
    my_keep_max=[ids_pass(maxima)',vals_smooth_pass(maxima)];
    
    if numel(my_keep_max)==0
        [maxima]=islocalmax(vals_smooth);
        new_max=[ids(maxima)',vals_smooth(maxima)];
        film_loc=new_max(new_max(:,2)==min(new_max(:,2)),1);
        mat_keep=mat_index_tot(film_loc(1,1),:);
    else        
        film_reduce=[(my_keep_max(1,1):my_keep_max(end,1))',vals_smooth(my_keep_max(1,1):my_keep_max(end,1))];
        [film_minima]=islocalmin(film_reduce(:,2));
        [film_maxima]=islocalmax(film_reduce(:,2));
        
        film_mins=film_reduce(film_minima,:);
        film_maxs=film_reduce(film_maxima,:);
        
        if and(numel(film_mins(:,2))==0,numel(film_maxs(:,2))==0)==1
            [film_maxima]=islocalmax(vals_smooth);
            val_keep=[ids(film_maxima)',vals_smooth(film_maxima)];
            ids_keep=val_keep(val_keep(:,2)==min(val_keep(:,2)),1);
            mat_keep=mat_index_tot(ids_keep,:);
        elseif and(numel(film_mins(:,2))==1,numel(film_maxs(:,2))==0)==1
            mat_keep=mat_index_tot(film_mins(1,1),:);
        elseif and(numel(film_mins(:,2))==0,numel(film_maxs(:,2))==1)==1
            mat_keep=mat_index_tot(film_maxs(1,1),:);
            %this means there is an air bubble in the middle of the film and the
            %grayscale values are such that there is a minimum where the middle of the
            %film should be
        elseif and(and(mod(numel(film_mins(:,2)),2)==1,mod(numel(film_maxs(:,2)),2)==0),numel(film_maxs(:,2))>0)==1
            minima_mid=film_mins(ceil(numel(film_mins(:,2))/2),:);
            mat_keep=mat_index_tot(minima_mid(1),:);
            %this is currently uncharted terriotry but I expect that there is only one
            %minima and no maxima.
        elseif and(and(mod(numel(film_mins(:,2)),2)==0,mod(numel(film_maxs(:,2)),2)==1),numel(film_mins(:,2))>0)==1
            maxima_mid=film_maxs(ceil(numel(film_maxs(:,2))/2),:);
            mat_keep=mat_index_tot(maxima_mid(1),:);
        else
            plot(ids,my_vals,'ob')
            hold on
            plot((1:numel(vals_smooth)),vals_smooth,'or')
            keyboard
        end
    end

    
    %Here we find the equation of the circle that fits the film. It is known
    %exactly if we know the height of the arc and two points on the arc. This
    %is known from th eperpendicular bisector we found earlier and from the
    %points of the vertices.
    ph_xy=mat_keep;
    %We want to find the center of the film that is either pinched up or pulled flat
    if del_x<0
        ph_xy(1)=floor(ph_xy(1));
    else
        ph_xy(1)=ceil(ph_xy(1));
    end
    if del_y<0
        ph_xy(2)=floor(ph_xy(2));
    else
        ph_xy(2)=ceil(ph_xy(2));
    end
    %we find the midpoint between the vertices and move it to account for
    %the error. Here we also account for moving the vertices in accordanc
    %with the error before we fit our circle.
    if errs==1
        p1_xy_err=[p1_x_pinch,p1_y_pinch];
        p2_xy_err=[p2_x_pinch,p2_y_pinch];
    else
        p1_xy_err=[p1_x_pull,p1_y_pull];
        p2_xy_err=[p2_x_pull,p2_y_pull];
    end
    %we need to know the distance between the two vertices
    d_verts=sqrt((p1_xy_err(1)-p2_xy_err(1))^2+(p1_xy_err(2)-p2_xy_err(2))^2);
    %We find the midpoint between the two vertices
    pm_xy=[(p1_xy_err(1)+p2_xy_err(1))/2,(p1_xy_err(2)+p2_xy_err(2))/2];
    h_circ=sqrt((ph_xy(1)-pm_xy(1)).^2+(ph_xy(2)-pm_xy(2)).^2);
    %We make sure we don't divide by zero
    if h_circ==0
        h_circ=1E-12;
    end
    rad_circ=d_verts^2/(8*h_circ)+h_circ;
    d_cen=rad_circ-h_circ;
    circ_plus_x=pm_xy(1)+d_cen/denom;
    circ_plus_y=pm_xy(2)+slope_bis*d_cen/denom;
    circ_minus_x=pm_xy(1)-d_cen/denom;
    circ_minus_y=pm_xy(2)-slope_bis*d_cen/denom;
    dist_plus_to_h=sqrt((circ_plus_x-ph_xy(1))^2+(circ_plus_y-ph_xy(2))^2);
    dist_minus_to_h=sqrt((circ_minus_x-ph_xy(1))^2+(circ_minus_y-ph_xy(2))^2);
    if dist_plus_to_h>dist_minus_to_h
        circ_cen_xy=[circ_plus_x,circ_plus_y];
    else
        circ_cen_xy=[circ_minus_x,circ_minus_y];
    end
    arc_eq_out(errs,:)=[circ_cen_xy,rad_circ];
    %Now we draw our circles. We have to test
    if xs_test>=ys_test
        xs_circ=(p1_xy_err(1):(p2_xy_err(1)-p1_xy_err(1))/(10*d_verts):p2_xy_err(1));
        if circ_cen_xy(2)>max([p1_xy_err(2),p2_xy_err(2)])
            ys_circ=circ_cen_xy(2)-sqrt(rad_circ^2-(xs_circ-circ_cen_xy(1)).^2);
        else
            ys_circ=circ_cen_xy(2)+sqrt(rad_circ^2-(xs_circ-circ_cen_xy(1)).^2);
        end
    else
        ys_circ=(p1_xy_err(2):(p2_xy_err(2)-p1_xy_err(2))/(10*d_verts):p2_xy_err(2));
        if circ_cen_xy(1)>max([p1_xy_err(1),p2_xy_err(1)])
            xs_circ=circ_cen_xy(1)-sqrt(rad_circ^2-(ys_circ-circ_cen_xy(2)).^2);
        else
            xs_circ=circ_cen_xy(1)+sqrt(rad_circ^2-(ys_circ-circ_cen_xy(2)).^2);
        end
    end
    %we once again build matrices depending on if we perform the pinch or
    %pull
    if errs==1
        xy_pinch=[xs_circ',ys_circ'];
        %The second time through the err loop we perform the "pull"
    else
        xy_pull=[xs_circ',ys_circ'];
    end
end


end