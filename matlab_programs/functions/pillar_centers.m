% +
% NAME: pillar_centers
%
% PURPOSE:%     
%
% CATEGORY:
%     stamping
%
% CALLING SEQUENCE:
%    vertex_locations=mc_vertex_2d(imiFileRead,imiFileWrite,walker_stack,kick_vec,thresh_vert)
%
% INPUTS: 
%    imiFileRead: Filename for the image we are trying to find the vertices
%    for.
%
%    imiFileWrite: Filename for the image with vertices displayed.
%
%    walker_stack: an NxNx120 matrix of veritex, where the third index
%    indicates the angle the primary leg makes with the x-axis.
%
%    kick_vec: size of allowed displacements and rotations for our Monte
%    Carlo method.
%
%    thresh_vert: is the maximum allowed energy for a walker to be  
%    considered a vertex.
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
%    written by: A. Chieco, UPenn, April 2021
%-
% 
function [vertex_loc,imi_out]=pillar_centers(struct_imi,walker_stack,centers,background)

imi_initial=struct_imi;
% create number of x,y,theta
num_x=numel(struct_imi(1,:));
num_y=numel(struct_imi(:,1));
num_theta=numel(walker_stack(1,1,:));
%The vertices are square so this is how big the square is
sz_vert=numel(walker_stack(:,1,1));

%we are making a list with centers, theta and number of pixels in the
%structuring element that are not background
xy_keep_tot=[centers,zeros(numel(centers(:,1)),2)];

%this is the total number of thetas we investigate
num_walkers=numel(xy_keep_tot(:,1));

for i1=1:num_walkers
    cen_inv=round(xy_keep_tot(i1,1:2));
    list_inv=zeros(num_theta,2);
    %this is the position in the image we are investigating
    my_pos=struct_imi(cen_inv(2)-floor(sz_vert/2):cen_inv(2)+floor(sz_vert/2),...
                      cen_inv(1)-floor(sz_vert/2):cen_inv(1)+floor(sz_vert/2));
    for i2=1:num_theta
        my_vertex=walker_stack(:,:,i2);  
        %We as a proxy for thr "energy" of the position we count the number
        %of non-background pixels that land within the structuring element
        mask_pos=my_pos.*my_vertex;
        %determine the energy of the theta position
        my_energy=numel(mask_pos(mask_pos<=background))-numel(my_vertex(my_vertex==0));
        list_inv(i2,:)=[i2,my_energy];
    end    
    list_keep=list_inv(list_inv(:,2)>0,:);
    if numel(list_keep)==0
        list_keep=list_inv(end,:);
    end
    if numel(list_keep(:,1))>1
        xy_keep_tot(i1,3:4)=list_keep(1,:);
    else
        xy_keep_tot(i1,3:4)=list_keep;
    end
end

vertex_loc=xy_keep_tot;

%now we make the final image of where our vertices are located
imi_out=struct_imi*0;

for i3=1:numel(vertex_loc(:,1))
    %We ensure our walkers do not overlap the boundary for x
    if vertex_loc(i3,1)>(num_x-floor(sz_vert/2))
        vertex_loc(i3,1)=num_x-floor(sz_vert/2);
    end
    if vertex_loc(i3,1)<floor(sz_vert/2)+1
        vertex_loc(i3,1)=floor(sz_vert/2)+1;
    end
    %We ensure our walkers do not overlap the boundary for y
    if vertex_loc(i3,2)>(num_y-floor(sz_vert/2))
        vertex_loc(i3,2)=num_y-floor(sz_vert/2);
    end
    if vertex_loc(i3,2)<floor(sz_vert/2)+1
        vertex_loc(i3,2)=floor(sz_vert/2)+1;
    end
    imi_square=imi_out(vertex_loc(i3,2)-floor(sz_vert/2):vertex_loc(i3,2)+floor(sz_vert/2),...
        vertex_loc(i3,1)-floor(sz_vert/2):vertex_loc(i3,1)+floor(sz_vert/2))+walker_stack(:,:,vertex_loc(i3,3));
    imi_out(vertex_loc(i3,2)-floor(sz_vert/2):vertex_loc(i3,2)+floor(sz_vert/2),...
        vertex_loc(i3,1)-floor(sz_vert/2):vertex_loc(i3,1)+floor(sz_vert/2))=imi_square;
end

imi_out(imi_out>0)=255;