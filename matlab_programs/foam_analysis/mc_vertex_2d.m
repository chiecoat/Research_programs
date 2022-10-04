% +
% NAME: mc_vertex_2d
%
% PURPOSE:%     
%
% CATEGORY:
%     Foam coarsening
%
% CALLING SEQUENCE:
%    vertex_locations=mc_vertex_2d(imiFileRead,imiFileWrite,walker_stack,kick_vec,thresh_vec)
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
%    thresh_vec: is a 2x1 vector. The first element is the minimum energy 
%    required for a "walker" to be placed onto the image, and the second 
%    element is the maximum allowed energy for a walker to be considered a 
%    vertex.
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
function vertex_loc=mc_vertex_2d(imiFileRead,FinalFileWrite,walker_stack,my_mask,kick_vec,thresh_vec)

foam_imi=imread(imiFileRead);
foam_imi=double(rgb2gray(foam_imi));
imi_initial=foam_imi;
imi_out=foam_imi;
num_x=numel(foam_imi(1,:));
num_y=numel(foam_imi(:,1));
num_theta=numel(walker_stack(1,1,:));
sz_vert=numel(walker_stack(:,1,1));


num_walkers=2E3;
%We start by seeding the locations of our vertices
vert_info_mat=zeros(num_walkers,5);
count_vec=(1:2E3);
my_count=1;
while_count=1;
while my_count<=num_walkers
   xy_seed=rand(3,1);
   x_y_theta=round([(num_x-sz_vert)*xy_seed(1),(num_y-sz_vert)*xy_seed(2),num_theta*xy_seed(3)]);
   x_y_theta=x_y_theta+1;
   if x_y_theta(3)>120
       x_y_theta(3)=1;
   end
   my_vertex=walker_stack(:,:,x_y_theta(3));
   my_pos=foam_imi(x_y_theta(2):x_y_theta(2)+sz_vert-1,x_y_theta(1):x_y_theta(1)+sz_vert-1);
   my_energy=sum(sum((my_pos.*my_mask-my_vertex).^2))/(sz_vert^2);
   if my_energy<=thresh_vec(1)
       vert_info_mat(count_vec(my_count),1:4)=[x_y_theta,my_energy];
       my_count=my_count+1;
       imi_initial(x_y_theta(2):x_y_theta(2)+sz_vert-1,x_y_theta(1):x_y_theta(1)+sz_vert-1)=my_vertex;
   end
   while_count=while_count+1;
   if while_count>1E5
       break
   end
end

%Now we perform our random kicks
for i1=1:1E3
    for i2=1:numel(vert_info_mat(:,1))
        xy_seed=rand(3,1);
        x_kick=2*kick_vec(1)*xy_seed(1)-kick_vec(1);
        y_kick=2*kick_vec(2)*xy_seed(2)-kick_vec(2);
        t_kick=2*kick_vec(3)*xy_seed(3)-kick_vec(3);
        new_pos=round([vert_info_mat(i2,1)+x_kick,vert_info_mat(i2,2)+y_kick,vert_info_mat(i2,3)+t_kick]);
        new_pos=new_pos+1;
        %We ensure our walkers do not overlap the boundary for x
        if new_pos(1)>(num_x-sz_vert)
            new_pos(1)=num_x-sz_vert;
        end
        if new_pos(1)<1
            new_pos(1)=1;
        end
        %We ensure our walkers do not overlap the boundary for y
        if new_pos(2)>(num_y-sz_vert)
            new_pos(2)=num_y-sz_vert;
        end
        if new_pos(2)<1
            new_pos(2)=1;
        end
        %we use the symmetry of the vertices to make sure everything is being
        %tracked between 0 and 119 degrese
        if new_pos(3)>120
            new_pos(3)=new_pos(3)-120;
        end
        if new_pos(3)<=0
            new_pos(3)=new_pos(3)+120;
        end
        my_vertex=walker_stack(:,:,new_pos(3));
        my_pos_new=foam_imi(new_pos(2):new_pos(2)+sz_vert-1,new_pos(1):new_pos(1)+sz_vert-1);
        my_energy_new=sum(sum((my_pos_new.*my_mask-my_vertex).^2))/(sz_vert^2);
        %This is the "zero temperature" condition
        if my_energy_new<=vert_info_mat(i2,4)
            vert_info_mat(i2,1:4)=[new_pos,my_energy_new];
        else
            vert_info_mat(5)=vert_info_mat(5)+1;
        end
    end 
end

vert_keep=vert_info_mat(vert_info_mat(:,4)<thresh_vec(2),:);
bark=0
for i3=1:numel(vert_keep(:,1))
     imi_out(vert_keep(i3,2):vert_keep(i3,2)+sz_vert-1,...
             vert_keep(i3,1):vert_keep(i3,1)+sz_vert-1)=walker_stack(:,:,vert_keep(i3,3));
end
imwrite(imi_out/255,FinalFileWrite);
% imwrite(imi_initial/255,InitialFileWrite)

vertex_loc=vert_keep(:,1:4);


end






















