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
%    None
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
%    written by: A. Chieco & Y. Shen, UPenn, September 2017
%-
function [arc_eq_out,x_circ_final,y_circ_final]=arc_find(p1_xy,p2_xy,foam_imi,bounds)

L_sys_x=numel(foam_imi(1,:));
L_sys_y=numel(foam_imi(:,1));


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
denom_factor=4;

%now we scan the image for points we want to try and fit circles to. The
%vertices are either "horizontally oriented" or "vertically oriented" and
%we have to distinguish between the two cases
xs_test=abs(p1_xy(1)-p2_xy(1));
ys_test=abs(p1_xy(2)-p2_xy(2));

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
bar=mean(my_vals);
val_keep=my_vals(my_vals<=bar); mat_keep=mat_index_tot(my_vals<=bar,:);
%This means that all values in val_keep and in my_vals are 0, so we just
%add one to all values of val_keep so we don't have to worry about
%weighting
if max(abs(val_keep-max(my_vals)))==0
    val_keep=val_keep+1;
end
ind_weighted=abs(val_keep-max(my_vals))/max(abs(val_keep-max(my_vals)));
x_weighted=sum(mat_keep(:,1).*ind_weighted)/sum(ind_weighted);
y_weighted=sum(mat_keep(:,2).*ind_weighted)/sum(ind_weighted);

% if and(and(x_weighted<=1479,x_weighted>=1414),and(y_weighted<=54,y_weighted>=15))==1
%     [x_weighted,y_weighted]
%     keyboard
% end

if or(or(p1_xy(1)==bounds(1,1),p1_xy(2)==bounds(2,1)),or(p2_xy(1)==bounds(1,2),p2_xy(2)==bounds(2,2)))==1
    xy_circ_fit=[p1_xy;p2_xy;x_weighted,y_weighted];
    arc_eq_out=CircleFitByPratt(xy_circ_fit);
    circ_cen_xy=arc_eq_out(1:2);
    rad_circ=arc_eq_out(3);
    %This means the sophisticated method did not work so we will use our
    %dumbed down method
    if isnan(rad_circ)==1
        ph_xy=[x_weighted,y_weighted];
        h_circ=sqrt((ph_xy(1)-pm_xy(1))^2+(ph_xy(2)-pm_xy(2))^2);
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
        arc_eq_out=[circ_cen_xy,rad_circ];
    end
else    
    ph_xy=[x_weighted,y_weighted];
    h_circ=sqrt((ph_xy(1)-pm_xy(1))^2+(ph_xy(2)-pm_xy(2))^2);
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
    arc_eq_out=[circ_cen_xy,rad_circ];
end
%Now we draw our circles on our images. We have to test the pixel
%values the circles encounter as they're being drawn.
if xs_test>=ys_test
    xs_keep=(p1_xy(1):(p2_xy(1)-p1_xy(1))/(10*d_verts):p2_xy(1))';
    if circ_cen_xy(2)>max([p1_xy(2),p2_xy(2)])
        ys_keep=circ_cen_xy(2)-sqrt(rad_circ^2-(xs_keep-circ_cen_xy(1)).^2);
    else
        ys_keep=circ_cen_xy(2)+sqrt(rad_circ^2-(xs_keep-circ_cen_xy(1)).^2);
    end
else
    ys_keep=(p1_xy(2):(p2_xy(2)-p1_xy(2))/(10*d_verts):p2_xy(2))';
    if circ_cen_xy(1)>max([p1_xy(1),p2_xy(1)])
        xs_keep=circ_cen_xy(1)-sqrt(rad_circ^2-(ys_keep-circ_cen_xy(2)).^2);
    else
        xs_keep=circ_cen_xy(1)+sqrt(rad_circ^2-(ys_keep-circ_cen_xy(2)).^2);
    end
end
% end

x_circ_final=xs_keep;
y_circ_final=ys_keep;


end
