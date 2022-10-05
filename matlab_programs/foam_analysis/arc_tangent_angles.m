% +
% NAME: arc_tangent_angles
%
% PURPOSE:%     
%
% CATEGORY:
%     Foam coarsening
%
% CALLING SEQUENCE:
%    tangents=arc_tangent_angles(verts,bubbles,arcs)
%
% INPUTS: 
%    verts = list of vertex locations and angles as well as error in the
%    locations and angles.verts will have the original vertex ID,
%    x,y,theta, location 'energy', bubble IDs for the original vertex (all will have 3)
%
%    bubbles = matrix with the bubble ID, the number of sides, the x
%    centroid, the y centroid, the vertex ids that make up the bubble.
%
%    arc_eqns = has the two vertex IDs, the x_center, the y_center and the
%    radius of the circle that makes up the film. 
%
%    ell = the distance from the point tangent to the circle that we want
%    our line segment to be. This number should not make a difference in
%    our calculations
%
% OPTIONAL INPUTS: (none)
%
% KEYWORD PARAMETERS: (none)
%
% OUTPUTS: a matrix whose elements are vertex id, neighbor id, 
%    (x,y) of the point tangent to the radius of the circle connecting the
%    vertex and its neighbor at the vertex, and the angle it makes with the
%    counterclockwise tangent line. The progrma also returns an image with
%    the tangent lines layed over the foam image.
% 
% SIDE EFFECTS: (none)
%
% MODIFICATION HISTORY:
%    written by: A. Chieco, UPenn, June 2018
%-
function [theta_tot,foam_arc_imi]=arc_tangent_angles(foam_arc_imi,verts,arcs,ell)

theta_tot=[-1,-1,-1,-1,-1];
for i1=1:numel(verts(:,1))
    vert_tan=verts(i1,:);
    v_id=vert_tan(1);
    xy_vert=vert_tan(2:3);
    v_check1=arcs(arcs(:,1)==v_id,:);
    v_check2=arcs(arcs(:,2)==v_id,:);
    v_check=[v_check1;v_check2];
    if numel(v_check(v_check(:,5)<0,5))>1
        continue
    end
    tangent_points=zeros(numel(v_check(:,1)),4);
    for i2=1:numel(v_check(:,1))
        if v_check(i2,5)<0
            tan_point=v_check(i2,3:4);
        else
            xy_circ=v_check(i2,3:4);
            rise=xy_circ(2)-xy_vert(2);
            if rise==0
                rise=1E-8;
            end
            run=xy_circ(1)-xy_vert(1);
            if abs(run)>1E8
                run=1E8*(run/abs(run));
            end
            tan_slope=-run/rise;
            y_int=xy_vert(2)-tan_slope*xy_vert(1);
            cap_A=1+tan_slope^2;
            cap_B=2*(tan_slope*y_int-xy_vert(1)-tan_slope*xy_vert(2));
            cap_C=y_int^2-2*xy_vert(2)*y_int+xy_vert(1)^2+xy_vert(2)^2-ell^2;
            root=cap_B^2-4*cap_A*cap_C;
            if root<0
                keyboard
                root=0;
            end
            x_plus=(-cap_B+sqrt(root))/(2*cap_A);
            y_plus=tan_slope*x_plus+y_int;
            x_minus=(-cap_B-sqrt(root))/(2*cap_A);
            y_minus=tan_slope*x_minus+y_int;
            n_id_tot=v_check(i2,1:2);
            n_id=n_id_tot(n_id_tot~=v_id);
            neighbor=verts(n_id,:);
            xy_neighbor=neighbor(2:3);
            d_vec=sqrt((xy_neighbor(1)-[x_plus,x_minus]).^2+(xy_neighbor(2)-[y_plus,y_minus]).^2);
            if d_vec(1)<d_vec(2)
                tan_point=[x_plus,y_plus];
            else
                tan_point=[x_minus,y_minus];
            end
        end
        tangent_points(i2,:)=[v_id,n_id,tan_point];
    end
    if numel(tangent_points(:,1))<=2
        point_order=[1,2];
    else       
        point_order=convhull(tangent_points(:,3),tangent_points(:,4));
    end
    theta_mat=zeros(numel(point_order)-1,5);
    if numel(tangent_points)==0
        continue
    end
    for i3=1:numel(point_order)-1        
        xy1=tangent_points(point_order(i3),3:4);
        vec1_scale=[xy1(1)-xy_vert(1),xy1(2)-xy_vert(2)];
        xy2=tangent_points(point_order(i3+1),3:4);
        vec2_scale=[xy2(1)-xy_vert(1),xy2(2)-xy_vert(2)];
        mag1=sqrt(sum(vec1_scale.^2));
        mag2=sqrt(sum(vec2_scale.^2));
        theta=acos(sum(vec1_scale.*vec2_scale)/(mag1*mag2));
        theta_mat(i3,:)=[tangent_points(point_order(i3+1),:),theta*180/pi];
    end
    theta_tot=[theta_tot;theta_mat];
end
theta_tot=theta_tot(2:end,:);

%We want to identify where there are bad vertex angles. This will put a
%mark over where there is a bad vertex
L_sys_y=numel(foam_arc_imi(:,1));
L_sys_x=numel(foam_arc_imi(1,:));

for i2=1:numel(theta_tot(:,1))
    vert_check=verts(theta_tot(i2,1),2:3);
    points=theta_tot(i2,3:4);
    x_diff=points(1)-vert_check(1);
    y_diff=points(2)-vert_check(2);
    slope=y_diff/x_diff;
    y_int=vert_check(2)-slope*vert_check(1);
    if abs(x_diff) > abs(y_diff)
        xs=floor((vert_check:x_diff/abs(x_diff):points(1)));
        ys=floor(slope*xs+y_int);        
    else
        ys=floor((vert_check(2):y_diff/abs(y_diff):points(2)));
        xs=floor((ys-y_int)/slope);
    end
    xy_mat=[xs',ys'];
    ys_keep=xy_mat(and(xy_mat(:,2)>=1,xy_mat(:,2)<=L_sys_y),:);
    xy_keep=ys_keep(and(ys_keep(:,1)>=1,ys_keep(:,1)<=L_sys_x),:);
    tan_ind=sub2ind([L_sys_y,L_sys_x],xy_keep(:,2),xy_keep(:,1));
    val=10*fix(abs(round(theta_tot(i2,5))-120)/10);
    foam_arc_imi(tan_ind)=val; 
end



