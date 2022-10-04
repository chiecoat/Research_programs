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
%    written by: A. Chieco, UPenn, October 2018
%-
clear all
close all


date='20181107';
savePath=['F:\Chieco\2d tracking\' date '-nitrogen-fill_height=79mm'];
imiPath=['F:\Chieco\2d tracking\' date '-nitrogen-fill_height=79mm\images'];
imi_list=(7176:7567);
for i1=1:numel(imi_list)
    im_num=imi_list(i1);%im_start+i1;
    im_num_str=num2str(im_num);
    im_name=['DSC_' im_num_str];
    %here we construct the matrix stack of vertices we will use to
    %reconstruct each foam image. We only need to do this once for th
    %ewhole stack
    if i1==1
        vert_name=im_name;
        walkerPath=['F:\foam coarsening\' date '-nitrogen-fill_height=79mm\DSC images\DSC_6017.jpg'];
        vertex_imi=imread(walkerPath);
        vertex_imi=double(rgb2gray(vertex_imi));
        %vertex_loc stores the x,y, and theta position of the vertex. Theta tracks
        %the angle of the right most leg of our vertex, which we will call the
        %primary leg
        box_size=211;
        radius=box_size/2;
        vertex_loc_clean=[3507,668,0];
        vertex_clean=vertex_imi(vertex_loc_clean(2):vertex_loc_clean(2)+box_size-1,vertex_loc_clean(1):vertex_loc_clean(1)+box_size-1);
        vertex_loc_bub=[2298,1750,105];
        vertex_bub=vertex_imi(vertex_loc_bub(2):vertex_loc_bub(2)+box_size-1,vertex_loc_bub(1):vertex_loc_bub(1)+box_size-1);
        vertex_loc_wet=[2963,2440,0];
        vertex_wet=vertex_imi(vertex_loc_wet(2):vertex_loc_wet(2)+box_size-1,vertex_loc_wet(1):vertex_loc_wet(1)+box_size-1);

        
        
        %Now we construct our circular mask
        x1 = box_size/2;
        y1 = box_size/2;
        v_rot=zeros(box_size,box_size,120,3);
        %Generate grid with binary mask representing the circle. 
        [xx,yy] = meshgrid((1:box_size)-y1,(1:box_size)-x1);
        mask=zeros(box_size,box_size);
        mask((xx.^2 + yy.^2)<radius^2)=1;
        for vert_type=1:3
            if vert_type==1
                vertex=vertex_clean;
            elseif vert_type==2
                vertex=vertex_bub;
            else
                vertex=vertex_wet;
            end
            %Now we fill in our matrix of rotated vertices
            t_start=-vertex_loc_clean(3); t_fin=120-vertex_loc_clean(3)-1;            
            for theta=t_start:t_fin
                v_new=imrotate(vertex,theta);
                vx_cen=numel(v_new(1,:))/2; vy_cen=numel(v_new(:,1))/2;
                v_new=v_new(vy_cen-floor(box_size/2):vy_cen+floor(box_size/2),...
                    vx_cen-floor(box_size/2):vx_cen+floor(box_size/2));
                v_rot(:,:,theta+vertex_loc_clean(3)+1,vert_type)=v_new.*mask;
            end
        end
    end
    vertexPath=['F:\foam coarsening\' date '-nitrogen-fill_height=79mm\skeleton images\' im_name ' skel.png' ];
    vertex_imi=imread(vertexPath);
    
    
    in=max([round(0.01*min(size(vertex_imi))),floor(box_size)]);
    if mod(in,2)>0
        in=in+1;
    end
    vertex_imi=vertex_imi(in/2:end-in/2,in/2:end-in/2);
    
    L_sys_y=numel(vertex_imi(:,1));
    L_sys_x=numel(vertex_imi(1,:));
    
    vertex_imi(vertex_imi<150)=0;
    vertex_imi(vertex_imi>150)=255;
    
    extra1=3;
    watershed_imi=zeros(L_sys_y+2*extra1,L_sys_x+2*extra1);
    %we are going to pad watershed imi and then draw a box around the foam
    watershed_imi(extra1:end-extra1-1,extra1:end-extra1-1)=vertex_imi;
    watershed_imi(1:extra1-1,1:end)=255;
    watershed_imi(end-extra1:end,1:end)=255;
    watershed_imi(1:end,1:extra1-1)=255;
    watershed_imi(1:end,end-extra1:end)=255;
    
    L_watershed_y=numel(watershed_imi(:,1));
    L_watershed_x=numel(watershed_imi(1,:));
    extra2=2*extra1;
    watershed_full=zeros(L_watershed_y+2*extra2,L_watershed_x+2*extra2);
    %we are going to pad watershed imi and then draw a box around the foam
    watershed_full(extra2:end-extra2-1,extra2:end-extra2-1)=watershed_imi;
    
    Ld_full=watershed(watershed_full);
    stats=regionprops(Ld_full,'Centroid');
    centroids=cat(1, stats.Centroid);
    
%     rgb = label2rgb(Ld_full,'jet',[.5 .5 .5]);
%     imshow(rgb);
    
    %we find the x and y of foam frame, plot the x and y in a figure
    xy_mat=[-1,-1];
    for i2=1:numel(Ld_full(1,:))
        ys=find(Ld_full(:,i2)==0);
        xy_new=zeros(numel(ys),2);
        xy_new(:,1)=i2;
        xy_new(:,2)=ys;
        xy_mat=[xy_mat;xy_new];
    end
    xy_mat=xy_mat(2:end,:);
    
    %here we find the location of the foam vertices
    xy_tag3=[xy_mat,zeros(numel(xy_mat(:,1)),3)];
    xy_tag3(:,5)=-1;
    for i2=1:numel(xy_mat(:,1))
        % establish a small box named Ld_small at each x and y. Because
        % of the padding of Ld_full this if statement is probably not used
        Ld_small=Ld_full((xy_mat(i2,2)-1):(xy_mat(i2,2)+1),(xy_mat(i2,1)-1):(xy_mat(i2,1)+1));
        % find all the different values in the small box
        unique_vals=unique(Ld_small);
        % find the non-zero value
        bub_id=find(unique_vals~=0);
        if and(numel(bub_id)==3,sum(unique_vals(bub_id)~=1)==3)==1
            % add the non-zero values into third through n columns of i2 row
            xy_tag3(i2,3:3+numel(bub_id)-1)=sort(unique_vals(bub_id));
        else
            continue
        end
    end
    % keep the first and second columns of xy_tag if the value of fifth column is greater than
    % zero
    xy_keep=xy_tag3(xy_tag3(:,5)>0,:);
    %We seed thevertex locatons as either the center
    xy_keep=[xy_keep(:,1:2);centroids(2:end,:)];

    
    %Now we give the vertics their new values so they
    xy_keep=[xy_keep(:,1)-extra2-extra1+1,xy_keep(:,2)-extra2-extra1+1];
    xy_keep(xy_keep(:,1)<=0,1)=1;
    xy_keep(xy_keep(:,2)<=0,2)=1;
    xy_keep(xy_keep(:,1)>=L_sys_x,1)=L_sys_x;
    xy_keep(xy_keep(:,2)>=L_sys_y,2)=L_sys_y;
    xy_keep_tot=[(1:numel(xy_keep(:,1)))',xy_keep(:,1)+in/2,xy_keep(:,2)+in/2,xy_keep(:,3:end)];
    clear xy_keep
    Ld=Ld_full(extra2+extra1:end-extra2-extra1-1,extra2+extra1:end-extra2-extra1-1);
%     rgb = label2rgb(Ld,'jet',[.5 .5 .5]);
%     imshow(rgb);
    
    foamPath=['F:\foam coarsening\' date '-nitrogen-fill_height=79mm\DSC images\' im_name '.jpg'];
    % read image file
    foam_pic=imread(foamPath);
    % convert image to gray image
    foam_pic=double(rgb2gray(foam_pic));
    %the threshol energy is going to be arbitrarily large the first
    %time we pass through the program. a good number is actually like
    %1.5E3
    thresh_vert=1.25E3;
    thresh_bub=9E3;
%     [vertex_locations,vertex_imi_final]=zoom_watershed_vertex_2d(foam_pic,v_rot,mask,xy_keep_tot,...
%                                         [box_size/20,box_size/20,120],[thresh_vert,thresh_bub],[box_size/10,box_size/2]);
% %     %vertex total will be a matrix that has the original vertex ID,
% %     %x,y,theta, location 'energy', bubble IDs for the original vertex (all will have 3)
%     vertex_total=[(1:numel(vertex_locations(:,1)))',vertex_locations(:,2:6)];
    vPath=[savePath '\vertices\vertices imi_' im_num_str];
    vertex_write=[vPath '_final.txt'];
%     dlmwrite(vertex_write,vertex_total,'delimiter',',','newline','pc','Precision','%1.8e')
%     vImPath=[savePath '\images\vertex final\vertex_overlay imi_' im_num_str];
%     imwrite(vertex_imi_final./255,[vImPath '_final.png']);
    vertex_total=dlmread(vertex_write);
%   %Here we find error in the locations of the vertices.    
    vert_with_errors=zoom_vertex_error(vertex_total,v_rot,mask,foam_pic,5,3);
    %this is just the uncertainty in x,y,theta
    vert_errors=vert_with_errors(:,end-2:end);
    %we have to connect the network of vertices. We have tracked one "leg"
    %of the vertices which will inform us which directions to look.
    vertex_total=zoom_vertex_connect(vertex_total,box_size);
     
    %Now we want to find all the vertices that belong to a bubble. We
    %already know the network of vertex connections. vertex total has
    %colums [id,x,y,theta,vert_type,energy,N_id_1,N_id_2,N_id_3] where N_id
    %is the neighbor ID of -1 if the vertex does not have 3 neighbors.
    %vertex_info will be returned with 
    %[id,x,y,theta,energy,bub_id_1,bub_id_2,bub_id_3,N_id_1,N_id_2,N_id_3]
    %and bubble info will have columns ...
    [bubble_info,vertex_info]=zoom_bubbles(vertex_total);
    bubPath=[savePath '\bubbles\bubbles imi_' im_num_str];
    bubbles_write=[bubPath '_final.txt'];
    dlmwrite(bubbles_write,bubble_info,'delimiter',',','newline','pc','Precision','%1.8e')
    %         plot(bubble_info(:,3),bubble_info(:,4),'ob')
    [arc_eqns,arc_errors,arc_imi]=zoom_arcs(foam_pic,vertex_info,vert_errors,box_size/2);
    arcPath=[savePath '\arcs\arcs imi_' im_num_str];
    arcs_write=[arcPath '_final.txt'];
    dlmwrite(arcs_write,arc_eqns,'delimiter',',','newline','pc','Precision','%1.8e')
    arcImPath=[savePath '\images\arcs overlay\arcs_overlay imi_' im_num_str];
    imwrite(arc_imi/255,[arcImPath '.png']);
    %Once we find the arcs of the bubbles we need to make sure that our
    %foam is well constructed. every vertex (except the ones on the
    %edges) need to be 3 fold connected and have an angles about 120
    %degrees. We are going to find vertices where this is not the case
    %and analyze them more critically to find all of the vertices in
    %that area.
    %             [tangents_and_thetas,foam_tan_imi]=arc_tangent_angles(foam_pic,vertex_info,arc_eqns,box_size/2);
    %             tanImPath=[savePath '\images\tangents overlay\tangents_overlay imi_' im_num_str];
    %             imwrite(foam_tan_imi/255,[tanImPath '.png']);
    %tangent lines is returned as a matrix whose elements are vertex
    %id, neighbor id, (x,y) of the point tangent to the radius of the
    %circle connecting the vertex and its neighbor at the vertex, and
    %the angle it makes with the counterclockwise tangent line.
    
    %We loop through the entire reconstruction again until we do not
    %have to update the vertex locations.
    %Once we are satisfied with our reconstruction
    my_bubble_stats=bubble_stats(vertex_total,bubble_info,arc_eqns,arc_errors);
    statsPath=[savePath '\stats\stats_bub_nxy_pace imi_' im_num_str];
    stats_write=[statsPath '_final.txt'];
    time_stamp=zeros(numel(my_bubble_stats(:,1)),1)+i1;
    dlmwrite(stats_write,[my_bubble_stats,time_stamp],'delimiter',',','newline','pc','Precision','%1.8e')
end

"all done"
