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
function [arc_eq_out,xy_pinch]=arc_error_fit(p1_xy,p2_xy,p1_err,p2_err,foam_imi)

L_sys_x=numel(foam_imi(1,:));
L_sys_y=numel(foam_imi(:,1));

%First we move the x positions of each vertex around. We either "pinch" or
%"pull" where we move the x positions closer together or further apart
if p1_xy(1)>p2_xy(1)
    p1_x_pinch=p1_xy(1)-p1_err(1);
    p2_x_pinch=p2_xy(1)+p2_err(1);
else
    p1_x_pinch=p1_xy(1)+p1_err(1);
    p2_x_pinch=p2_xy(1)-p2_err(1);
end
%We do the same thing for the y positions
if p1_xy(2)>p2_xy(2)
    p1_y_pinch=p1_xy(2)-p1_err(2);
    p2_y_pinch=p2_xy(2)+p2_err(2);
else
    p1_y_pinch=p1_xy(2)+p1_err(2);
    p2_y_pinch=p2_xy(2)-p2_err(2);
end

del_x=abs(p1_err(1)-p2_err(1));
del_y=abs(p1_err(2)-p2_err(2));
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
denom_factor=4;

%"Horizontally oriented" vertices lie further apart in x than y
if xs_test>=ys_test
    %These two points make a line that bisects the the vertices
    x_seek_plus=pm_xy(1)+d_verts/(denom_factor*denom);
    x_seek_minus=pm_xy(1)-d_verts/(denom_factor*denom);
    %There is a chance the two points are colinear
    if x_seek_plus~=x_seek_minus
        xs_seek=(x_seek_plus:(x_seek_minus-x_seek_plus)/d_verts:x_seek_minus);
        ys_seek=slope_bis*xs_seek+y_int_bis;
    else
        ys_seek=(-d_verts/denom_factor:1:d_verts/denom_factor)+pm_xy(2);
        xs_seek=zeros(1,numel(ys_seek))+x_seek_plus;
    end
else
%"Vertically oriented" vertices lie further apart in y than x
    %These two points make a line that bisects the the vertices
    y_seek_plus=pm_xy(2)+slope_bis*d_verts/(denom_factor*denom);
    y_seek_minus=pm_xy(2)-slope_bis*d_verts/(denom_factor*denom);
    %There is a chance the two points are colinear
    if y_seek_plus~=y_seek_minus
        ys_seek=(y_seek_plus:(y_seek_minus-y_seek_plus)/d_verts:y_seek_minus);
        xs_seek=(ys_seek-y_int_bis)/slope_bis;
    else
        xs_seek=(-d_verts/denom_factor:1:d_verts/denom_factor)+pm_xy(1);
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
% valley_loc_tot=find(my_vals==min(my_vals));
% valley_loc=valley_loc_tot(1);
% film_cen=val_ind(valley_loc);
% [ph_y,ph_x]=ind2sub([L_sys_y,L_sys_x],film_cen);
bar=mean(my_vals);
val_keep=my_vals(my_vals<=bar); mat_keep=mat_index_tot(my_vals<=bar,:);
%This means that all values in val_keep and in my_vals are 0, so we just
%add one to all values of val_keep so we don't have to worry about
%weighting
if max(abs(val_keep-max(my_vals)))==0
    val_keep=val_keep+1;
end
ind_weighted=abs(val_keep-max(my_vals))/max(abs(val_keep-max(my_vals)));
ph_x=sum(mat_keep(:,1).*ind_weighted)/sum(ind_weighted);
ph_y=sum(mat_keep(:,2).*ind_weighted)/sum(ind_weighted);

%we find the midpoint between the vertices and move it to account for
%the error. Here we also account for moving the vertices in accordanc
%with the error before we fit our circle.
p1_xy_err=[p1_x_pinch,p1_y_pinch];
p2_xy_err=[p2_x_pinch,p2_y_pinch];
%Now we redo some of the other parameters we need
slope_m=(p2_xy_err(2)-p1_xy_err(2))/(p2_xy_err(1)-p1_xy_err(1));
%Slope 0 will give us an infinite radius, so we consider vertices that lie
%on in the same plane
if slope_m==0
    slope_m=1E-12;
end
%We find the line that perpendicularly bisects the two vertices
slope_bis_err=-1/slope_m;
denom_err=sqrt(1+slope_bis_err^2);
%we need to know the distance between the two vertices
d_verts=sqrt((p1_xy_err(1)-p2_xy_err(1))^2+(p1_xy_err(2)-p2_xy_err(2))^2);
%We find the midpoint between the two vertices
pm_xy=[(p1_xy_err(1)+p2_xy_err(1))/2,(p1_xy_err(2)+p2_xy_err(2))/2];
if ph_x>pm_xy(1)
    x_sign=1;
else
    x_sign=-1;
end
if ph_y>pm_xy(2)
    y_sign=1;
else
    y_sign=-1;
end
ph_xy=[ph_x(1)+x_sign*del_x,ph_y(1)+y_sign*del_y];
%We want to find the center of the film that is either pinched up or pulled flat
% if x_sign*del_x<0
%     ph_xy(1)=floor(ph_xy(1));
% else
%     ph_xy(1)=ceil(ph_xy(1));
% end
% if y_sign*del_y<0
%     ph_xy(2)=floor(ph_xy(2));
% else
%     ph_xy(2)=ceil(ph_xy(2));
% end
h_circ=sqrt((ph_xy(1)-pm_xy(1)).^2+(ph_xy(2)-pm_xy(2)).^2);
%We make sure we don't divide by zero
if h_circ==0
    h_circ=1E-12;
end
rad_circ=d_verts^2/(8*h_circ)+h_circ;
d_cen=rad_circ-h_circ;
circ_plus_x=pm_xy(1)+d_cen/denom_err;
circ_plus_y=pm_xy(2)+slope_bis_err*d_cen/denom_err;
circ_minus_x=pm_xy(1)-d_cen/denom_err;
circ_minus_y=pm_xy(2)-slope_bis_err*d_cen/denom_err;
dist_plus_to_h=sqrt((circ_plus_x-ph_xy(1))^2+(circ_plus_y-ph_xy(2))^2);
dist_minus_to_h=sqrt((circ_minus_x-ph_xy(1))^2+(circ_minus_y-ph_xy(2))^2);
if dist_plus_to_h>dist_minus_to_h
    circ_cen_xy=[circ_plus_x,circ_plus_y];
else
    circ_cen_xy=[circ_minus_x,circ_minus_y];
end
arc_eq_out=[circ_cen_xy,rad_circ];
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
xy_pinch=[xs_circ',ys_circ'];




end