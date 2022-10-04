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
function [arc_eq_out,x_circ_final,y_circ_final]=zoom_arc_find(p1_xy,p2_xy,foam_imi,bub_thresh,film_thresh)

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
denom_factor=2;

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
ids=(1:numel(my_vals));
vals_first_pass=my_vals(my_vals<=bub_thresh); 
ids_first_pass=ids(my_vals<=bub_thresh);
%now we identify the region between the films that we want to fit a circle
%to. We find the location of the local maxima 


plot((1:numel(vals_first_pass)),vals_first_pass,'or')
if and(and(p1_xy(1)==2058,p1_xy(2)==560),and(p2_xy(1)==2336,p2_xy(2)==454))==1
    keyboard
end

[minima]=islocalmin(vals_first_pass);
%We only keep minima that fall below a threshold we decide determined where
%a film is and it is not just noise in the data
vals_second_pass=vals_first_pass(minima);
ids_second_pass=ids_first_pass(minima);
bounds=find(vals_second_pass<film_thresh);
valley_min=ids_second_pass(bounds(1));
valley_max=ids_second_pass(bounds(end));

%these are the pixels that land within our film
val_width=my_vals(valley_min:valley_max);
ids_width=(valley_min:valley_max);

ids_in=ids_width(val_width>film_thresh);
id_diff=abs(ids_in(1:end-1)-ids_in(2:end));
id_edge_one=find(id_diff==max(id_diff));
if numel(id_edge_one)==0
    ids_in=ids_width;
    id_edge=[1,numel(ids_width)];
elseif numel(id_edge_one)==1
    id_diff(id_edge_one)=0;
    id_edge_two=find(id_diff==max(id_diff));
    id_edge=sort([id_edge_one,id_edge_two]);
    id_edge=[id_edge(1)+1,id_edge(2)];    
elseif numel(id_edge_one)==2
    id_edge=[id_edge_one(1)+1,id_edge_one(2)];    
elseif numel(id_edge_one)>2
    id_edge=[id_edge_one(1)+1,id_edge_one(end)];
end

if id_edge(1)~=id_edge(2)
    val_keep=my_vals(ids_in(id_edge(1)):ids_in(id_edge(2)));
    mat_keep=mat_index_tot(ids_in(id_edge(1)):ids_in(id_edge(2)),:);
else
    val_keep=my_vals(ids_width(1):ids_width(end));
    mat_keep=mat_index_tot(ids_width(1):ids_width(end),:);  
end
    

%Now we estimate the location we need to fit our cirlce to.
ind_weighted=abs(val_keep-min(my_vals))/max(abs(val_keep-min(my_vals)));
x_weighted=sum(mat_keep(:,1).*ind_weighted)/sum(ind_weighted);
y_weighted=sum(mat_keep(:,2).*ind_weighted)/sum(ind_weighted);

%Here we find the equation of the circle that fits the film. It is known
%exactly if we know the height of the arc and two points on the arc. This
%is known from th eperpendicular bisector we found earlier and from the
%points of the vertices.
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