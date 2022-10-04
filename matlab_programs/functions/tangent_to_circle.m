%+
%NAME: tangent_to_circle
%
%PURPOSE: 
%    Surface Evolver creates a foam where the vertices and plateau borders
%    obey periodic boundary conditions. We have adjusted the vertices 
%
%CATEGORY:
%    Foam Coarsening: Simulations
%
%CALLING SEQUENCE:
%   
%INPUTS:        
% 
%OPTIONAL INPUTS: None
%
%KEYWORD PARAMETERS: None
% 
%OUTPUTS: 
% 
%  SIDE EFFECTS: (none)
% 
%  MODIFICATION HISTORY:
%     written by: A. Chieco, UPenn, May 2016
%--------------------------------------------------------------------------
function tangent_return=tangent_to_circle(xy_in,circ_cen)

%this is hte vector pointing from the center of the circle to the point
%we want to find a tangent line for
RV_x=xy_in(1,1)-circ_cen(1);
RV_y=xy_in(1,2)-circ_cen(2);

tan_mag=sqrt((xy_in(1,1)-xy_in(2,1))^2+(xy_in(1,2)-xy_in(2,2))^2);

RV_mag=sqrt(RV_x^2+RV_y^2);
denom=RV_x/RV_y+RV_y/RV_x;

tan_x=-RV_mag*tan_mag/(RV_x*denom);
tan_y=RV_mag*tan_mag/(RV_y*denom);

%THere are two tangent lines at this point on the circle. They are
%separated by a 180 degree rotation. Here we find the other tangent line
tan_x_rot=-tan_x;
tan_y_rot=-tan_y;
           
tan_pos=[tan_x+xy_in(1,1),tan_y+xy_in(1,2); tan_x_rot+xy_in(1,1),tan_y_rot+xy_in(1,2)]; 
            
dist_mat=[sqrt((tan_pos(1,1)-xy_in(2,1))^2+(tan_pos(1,2)-xy_in(2,2))^2), ...
          sqrt((tan_pos(2,1)-xy_in(2,1))^2+(tan_pos(2,2)-xy_in(2,2))^2)]; 
     
tangent_return=tan_pos(find(dist_mat==min(dist_mat)),:);  
            
          
           
end