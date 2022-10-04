% +
% NAME: mc_vertex_2d
%
% PURPOSE:%     
%
% CATEGORY:
%     Foam coarsening
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
%    written by: A. Chieco & Y. Shen, UPenn, September 2017
%-
% 
function [vertex_loc,imi_out]=zoom_watershed_vertex_2d(foam_imi,walker_stack,my_mask,xy_tot,kick_vec,thresh_vec,dist_thresh)

imi_initial=foam_imi;
% create number of x,y,theta
num_x=numel(foam_imi(1,:));
num_y=numel(foam_imi(:,1));
num_theta=numel(walker_stack(1,1,:,1));
%The vertices are square so this is how big the square is
sz_vert=numel(walker_stack(:,1,1,1));

x_bound_min=floor(sz_vert/2)+1;
x_bound_max=num_x-floor(sz_vert/2);
y_bound_min=floor(sz_vert/2)+1;
y_bound_max=num_y-floor(sz_vert/2);

xy_keep_tot=xy_tot(:,2:3);

num_walkers=numel(xy_keep_tot(:,1));
%We start by seeding the locations of our vertices
vert_info_mat=zeros(num_walkers,6);

for my_count=1:num_walkers
    %xy_seed is just one row at each time,first column is x of vertice,
    %second column is y of vertice, third column is a random number between
    %0 and 1. We assume all vertices look good.
    vert_type=1;
    xy_seed=[xy_keep_tot(my_count,1:2),rand(1,1)];
    %make the elements of xy_seed be the closest integer
    x_y_theta_type=[round([xy_seed(1),xy_seed(2),num_theta*xy_seed(3)+1]),vert_type];
    %theta are symmetric around 120 degrees. We only need theta numbers that
    %are between 0-119. If theta=120, the images is the same as theta=0, so
    %we rest the condition here.
    if x_y_theta_type(3)>120
        x_y_theta_type(3)=1;
    end
    my_vertex=walker_stack(:,:,x_y_theta_type(3),1);
    %we make sure the vertex stays within the image, for the x boundaries
    if x_y_theta_type(1)<=x_bound_min
        x_y_theta_type(1)=x_bound_min;
    end
    if x_y_theta_type(1)>=x_bound_max
        x_y_theta_type(1)=x_bound_max;
    end
    %we check again for the y boundaries
    if x_y_theta_type(2)<=y_bound_min
        x_y_theta_type(2)=y_bound_min;
    end
    if x_y_theta_type(2)>=y_bound_max
        x_y_theta_type(2)=y_bound_max;
    end
    my_pos=foam_imi(x_y_theta_type(2)-floor(sz_vert/2):x_y_theta_type(2)+floor(sz_vert/2),...
        x_y_theta_type(1)-floor(sz_vert/2):x_y_theta_type(1)+floor(sz_vert/2));
    % We need to sum the values in the "my_pos" image to get an energy value
    my_energy=sum(sum((my_pos.*my_mask-my_vertex).^2))/(sz_vert^2);
    % determine the initial of vertex
    vert_info_mat(my_count,1:5)=[x_y_theta_type,my_energy];
    imi_initial(x_y_theta_type(2)-floor(sz_vert/2):x_y_theta_type(2)+floor(sz_vert/2),...
        x_y_theta_type(1)-floor(sz_vert/2):x_y_theta_type(1)+floor(sz_vert/2))=my_vertex;
end
remove=find(vert_info_mat(:,5)<thresh_vec(2));
vert_info_mat=vert_info_mat(remove,:);
xy_keep_tot=xy_keep_tot(remove,:);
xy_calc=xy_tot(remove,2:3);

%each vertex is allowed to explore the energy landscape that is as big as
%half the distance to its neighboring vertex. This way vertices do not
%accidentally swap places
dist_thresh_vec=zeros(numel(xy_keep_tot(:,1)),1);
for dist_calc=1:numel(xy_keep_tot(:,1))
    dist_vec=sort(sqrt((xy_keep_tot(dist_calc,1)-xy_calc(:,1)).^2+(xy_keep_tot(dist_calc,2)-xy_calc(:,2)).^2));
    %The first value of dist_vec is always zero
    dist_thresh_vec(dist_calc)=dist_vec(2)/3;
    %this is only possible if two vertices overlap at an edge. This is
    %purly coincidental but necessary to differentiate between artificial
    %vertices at an edge and actual vertices on a boudary between bubbles.
    %So we keep both but we kow the artifical vertex because it overlapy
    %with the cell_ID=1. That one stays at the edge and the other vertex
    %gets bumped off but only byb one pixel
        if numel(dist_vec(dist_vec==0))>1
            %we check if we are on an x bounadry
            if xy_keep_tot(dist_calc,1)<=x_bound_min
                xy_keep_tot(dist_calc,1)=x_bound_min+1;
            end
            if xy_keep_tot(dist_calc,1)>=x_bound_max
                xy_keep_tot(dist_calc,1)=x_bound_min-1;
            end
            %we check again for the y boundaries
            if xy_keep_tot(dist_calc,2)<=y_bound_min
                xy_keep_tot(dist_calc,2)=y_bound_min+1;
            end
            if xy_keep_tot(dist_calc,2)>=y_bound_max
                xy_keep_tot(dist_calc,2)=y_bound_min-1;
            end
        end
end
dist_thresh_vec(dist_thresh_vec>dist_thresh(2))=dist_thresh(2);
dist_thresh_vec(dist_thresh_vec<dist_thresh(1))=dist_thresh(1);

%Now we perform our random kicks
for i1=1:5E3
%     imi_mc_step=foam_imi;
    xy_seed=rand(numel(vert_info_mat(:,1)),3);
    for i2=1:numel(vert_info_mat(:,1))
        %we generate random numbers for the random kicks. The size of the
        %kicks are set by the elements in kick_vec (which is an input)
        x_kick=round(2*kick_vec(1)*xy_seed(i2,1)-kick_vec(1));
        y_kick=round(2*kick_vec(2)*xy_seed(i2,2)-kick_vec(2));
        t_kick=round(2*kick_vec(3)*xy_seed(i2,3)-kick_vec(3));
        new_type=floor(3*rand(1,1))+1;
        if new_type>3
            new_type=3;
        end
        new_pos=round([vert_info_mat(i2,1)+x_kick,vert_info_mat(i2,2)+y_kick,vert_info_mat(i2,3)+t_kick]);
        %We ensure our walkers do not overlap the boundary for x
        if new_pos(1)<=x_bound_min
            new_pos(1)=x_bound_min;
        end
        if new_pos(1)>=x_bound_max
            new_pos(1)=x_bound_max;
        end
        %We ensure our walkers do not overlap the boundary for y
        if new_pos(2)<=y_bound_min
            new_pos(2)=y_bound_min;
        end
        if new_pos(2)>=y_bound_max
            new_pos(2)=y_bound_max;
        end
        %we use the symmetry of the vertices to make sure everything is being
        %tracked between 0 and 119 degree
        if new_pos(3)>120
            new_pos(3)=new_pos(3)-120;
        end
        if new_pos(3)<=0
            new_pos(3)=new_pos(3)+120;
        end
        dist_from_start=sqrt((new_pos(1)-xy_keep_tot(i2,1))^2+(new_pos(2)-xy_keep_tot(i2,2))^2);
        if dist_from_start>dist_thresh_vec(i2)
            my_energy_new=1E5;
        else
            if or(new_pos(3)>120,new_pos(3)<=0)
                keyboard
            end
            my_vertex=walker_stack(:,:,new_pos(3),new_type);
            my_pos_new=foam_imi(new_pos(2)-floor(sz_vert/2):new_pos(2)+floor(sz_vert/2),...
                new_pos(1)-floor(sz_vert/2):new_pos(1)+floor(sz_vert/2));
            my_energy_new=sum(sum((my_pos_new.*my_mask-my_vertex).^2))/(sz_vert^2);
        end
        %The first statement in the or is a "zero temperature" condition
        if my_energy_new<=vert_info_mat(i2,5)
            vert_info_mat(i2,1:5)=[new_pos,new_type,my_energy_new];
        else
            vert_info_mat(i2,6)=vert_info_mat(i2,6)+1;
        end
%         %We add some temperature to kick vertices out of local minima.
%         if and(vert_info_mat(i2,5)>=thresh_vec(1),my_energy_new<=2*thresh_vec(1))      
%             vert_info_mat(i2,1:5)=[new_pos,new_type,my_energy_new];
%         end
%         if or(i1==1,mod(i1,40)==0)==1
%             imi_mc_step(vert_info_mat(i2,2)-floor(sz_vert/2):vert_info_mat(i2,2)+floor(sz_vert/2),...
%             vert_info_mat(i2,1)-floor(sz_vert/2):vert_info_mat(i2,1)+floor(sz_vert/2))=walker_stack(:,:,vert_info_mat(i2,3));
%         end
    end
    %We can write out each stpe to show a movie of the monte carlo steps
    if or(i1==1,mod(i1,40)==0)==1
        i1
%         imwrite(imi_mc_step/255,[FileWrite '_mcstep_' num2str(i1) '.png']);
    end
end
i1

%we add back the IDs for the vertices
vert_keep1=vert_info_mat(vert_info_mat(:,5)<thresh_vec(1),1:5);
vert_keep1=[(1:numel(vert_keep1(:,1)))',vert_keep1];
% remove the superfluous vertices. The return matrix will include a vector
% of ids in the first column. So vert_keep has the form
%[ids,x,y,theta,E]
vert_keep=remove_duplicates([vert_keep1(:,1:4),vert_keep1(:,6)],sz_vert/4);

vertex_loc=vert_keep1(vert_keep(:,1),:);

%now we make the final image of where our vertices are located
imi_out=foam_imi;
for i3=1:numel(vertex_loc(:,1))
    %We ensure our walkers do not overlap the boundary for x
    if vertex_loc(i3,2)>(num_x-floor(sz_vert/2))
        vertex_loc(i3,2)=num_x-floor(sz_vert/2);
    end
    if vertex_loc(i3,2)<floor(sz_vert/2)+1
        vertex_loc(i3,2)=floor(sz_vert/2)+1;
    end
    %We ensure our walkers do not overlap the boundary for y
    if vertex_loc(i3,3)>(num_y-floor(sz_vert/2))
        vertex_loc(i3,3)=num_y-floor(sz_vert/2);
    end
    if vertex_loc(i3,3)<floor(sz_vert/2)+1
        vertex_loc(i3,3)=floor(sz_vert/2)+1;
    end
    %we check what theta our vertices are at as well
    if vertex_loc(i3,4)>120
        vertex_loc(i3,4)=vertex_loc(i3,4)-120;
    end
    if vertex_loc(i3,4)<=0
        vertex_loc(i3,4)=vertex_loc(i3,4)+120;
    end
    imi_out(vertex_loc(i3,3)-floor(sz_vert/2):vertex_loc(i3,3)+floor(sz_vert/2),...
        vertex_loc(i3,2)-floor(sz_vert/2):vertex_loc(i3,2)+floor(sz_vert/2))=walker_stack(:,:,vertex_loc(i3,4),vertex_loc(i3,5));
end



    
    
    

