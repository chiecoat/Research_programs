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
%    written by: A. Chieco, UPenn, September 2017
%-
clear all 
close all


for folders=2:3
%     if folders==1
%         savePath='F:\Chieco\2d tracking\20170906-nitrogen-fill_height=48mm';
%         imiPath=[savePath '\images'];        
%         coarsenPath='F:\foam coarsening\20170906-nitrogen-fill_height=48mm';
%         folder_path='106NCD80';
%         im_start=400;
%         i_start=67;
%         i_fin=770;
%         step=1;
%     end
    if folders==2
        savePath='F:\Chieco\2d tracking\20170907-nitrogen-fill_height=48mm';
        imiPath=[savePath '\images'];        
        coarsenPath='F:\foam coarsening\20170907-nitrogen-fill_height=48mm';
        folder_path='107NCD80';
        im_start=400;
        i_start=591;
        i_fin=706;
        step=1;        
    end
    if folders==3
        savePath='F:\Chieco\2d tracking\20170908-nitrogen-fill_height=48mm';
        imiPath=[savePath '\images'];        
        coarsenPath='F:\foam coarsening\20170908-nitrogen-fill_height=48mm';
        folder_path='103NCD80';
        im_start=400;
        i_start=0;
        i_fin=678;
        step=1;        
    end        
    for i1=i_start:step:i_fin
        im_num=im_start+i1;
        im_num_str=num2str(im_num);
        if im_num<1E3
            im_name=['DSC_0' im_num_str];
        else
            im_name=['DSC_' im_num_str];
        end
        %here we construct the matrix stack of vertices we will use to
        %reconstruct each foam image. We only need to do this once for th
        %ewhole stack
        if i1==i_start
            vert_name=im_name;
            walkerPath=['F:\foam coarsening\20170906-nitrogen-fill_height=48mm\106NCD80\DSC_0600.JPG'];
            vertex_imi=imread(walkerPath);
            vertex_imi=double(rgb2gray(vertex_imi));
            %vertex_loc stores the x,y, and theta position of the vertex. Theta tracks
            %the angle of the right most leg of our vertex, which we will call the
            %primary leg
            vertex_loc=[1131,863,0];
            box_size=23;
            vertex=vertex_imi(vertex_loc(2):vertex_loc(2)+box_size-1,vertex_loc(1):vertex_loc(1)+box_size-1);
            %Now we construct our circular mask
            x1 = box_size/2;
            y1 = box_size/2;
            radius = box_size/2;
            %Generate grid with binary mask representing the circle. Credit to Jonas for original code.
            [xx,yy] = meshgrid((1:box_size)-y1,(1:box_size)-x1);
            mask=zeros(box_size,box_size);
            mask((xx.^2 + yy.^2)<radius^2)=1;
            %Now we fill in our matrix of rotated vertices
            t_start=-vertex_loc(3); t_fin=120-vertex_loc(3)-1;
            v_rot=zeros(box_size,box_size,120);
            for theta=t_start:t_fin
                v_new=imrotate(vertex,theta);
                vx_cen=numel(v_new(1,:))/2; vy_cen=numel(v_new(:,1))/2;
                v_new=v_new(vy_cen-box_size/2+1:vy_cen+box_size/2,vx_cen-box_size/2+1:vx_cen+box_size/2);
                v_rot(:,:,theta+vertex_loc(3)+1)=v_new.*mask;
            end
        end
        vertexPath=[coarsenPath '\skeleton images\skeleton ' im_name '.png' ];
        vertex_imi=imread(vertexPath);
        L_sys_x=numel(vertex_imi(1,:));
        L_sys_y=numel(vertex_imi(:,1));
        
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
        
        L_watershed_x=numel(watershed_imi(1,:));
        L_watershed_y=numel(watershed_imi(:,1));
        extra2=2*extra1;
        watershed_full=zeros(L_watershed_y+2*extra2,L_watershed_x+2*extra2);
        %we are going to pad watershed imi and then draw a box around the foam
        watershed_full(extra2:end-extra2-1,extra2:end-extra2-1)=watershed_imi;
        
        Ld_full=watershed(watershed_full);
        Ld=Ld_full(extra2+extra1:end-extra2-extra1-1,extra2+extra1:end-extra2-extra1-1);
        rgb = label2rgb(Ld,'jet',[.5 .5 .5]);
        %     imshow(rgb);
        % we find the x and y of foam frame, plot the x and y in a figure
        thresh = find(Ld==0);
        xy_mat=[-1,-1];
        for i2=1:numel(Ld(1,:));
            ys=find(Ld(:,i2)==0);
            xy_new=zeros(numel(ys),2);
            xy_new(:,1)=i2;
            xy_new(:,2)=ys;
            xy_mat=[xy_mat;xy_new];
        end
        xy_mat=xy_mat(2:end,:);
%         plot(xy_mat(:,1),xy_mat(:,2),'.r');
%         hold on
        
        
        xy_tag=[xy_mat,zeros(numel(xy_mat(:,1)),3)];
        xy_tag(:,5)=-1;
        for i3=1:numel(xy_mat(:,1))
            % establish a small box named Ld_small at each x and y
            if or(xy_mat(i3,2)-1<1,xy_mat(i3,2)+1>numel(Ld(:,1)))==1
                continue
            elseif or(xy_mat(i3,1)-2<1,xy_mat(i3,1)+1>numel(Ld(1,:)))==1
                continue
            end
            Ld_small=Ld((xy_mat(i3,2)-1):(xy_mat(i3,2)+1),(xy_mat(i3,1)-1):(xy_mat(i3,1)+1));
            % find all the different values in the small box
            unique_vals=unique(Ld_small);
            % find the non-zero value
            bub_id=find(unique_vals~=0);
            % add the non-zero values into third to fifth columns of i2 row
            if numel(bub_id)>3
                continue
            end
            xy_tag(i3,3:3+numel(bub_id)-1)=unique_vals(bub_id);
        end
        % keep the first and second columns of xy_tag if the value of fifth column is greater than
        % zero
        xy_keep=xy_tag(xy_tag(:,5)>0,:);
        %We can adjust "xy_keep_tot" in order to add other locations along the edge
        %of the foam, right now it only considers vertices in the bulk
        xy_keep_tot=[(1:numel(xy_keep(:,1)))',xy_keep];
        
        foamPath=[coarsenPath '\' folder_path '\' im_name '.JPG'];
        vImPath=[imiPath '\vertex final\vertex_overlay imi_' im_num_str];
        vPath=[savePath '\vertices\vertices imi_' im_num_str];
        % read image file
        foam_pic=imread(foamPath);
        % convert image to gray image
        foam_pic=double(rgb2gray(foam_pic));
        vertex_locations=watershed_vertex_2d(foam_pic,vImPath,v_rot,mask,xy_keep_tot,[1,1,30],1.5E3);
%         %vertex total will be a matrix that has the original vertex ID,
%         %x,y,theta, location 'energy', bubble IDs for the original vertex (all will have 3)
        vertex_total=[(1:numel(vertex_locations(:,1)))',vertex_locations(:,2:5),xy_keep_tot(vertex_locations(:,1),4:6)];
        vertex_write=[vPath '_final.txt'];
        dlmwrite(vertex_write,vertex_total,'delimiter',',','newline','pc','Precision','%1.8e')
%         vertex_total=dlmread(vertex_write);
        %Here we find error in the locations of the vertices.
        vert_with_errors=vertex_error(vertex_total,v_rot,mask,foam_pic,5,3);
        %this is just the uncertainty in x,y,theta
        vert_errors=vert_with_errors(:,end-2:end);
        %Now we want to find all the vertices that belong to a bubble. We
        %already know which bubble each vertex belong to because they were
        %tagged in the watershed
        %     hold on
        %     plot(vertex_total(:,2),vertex_total(:,3),'ob')
        %     hold off
        bubble_info=watershed_bubbles(vertex_total);
        bubPath=[savePath '\bubbles\bubbles imi_' im_num_str];
        bubbles_write=[bubPath '_final.txt'];
        dlmwrite(bubbles_write,bubble_info,'delimiter',',','newline','pc','Precision','%1.8e')
%         plot(bubble_info(:,3),bubble_info(:,4),'ok')
        [arc_eqns,arc_errors,arc_imi]=watershed_arcs(foam_pic,bubble_info,vertex_total,vert_errors);
        arcPath=[savePath '\arcs\arcs imi_' im_num_str];
        arcs_write=[arcPath '_final.txt'];
        arcImPath=[imiPath '\arcs overlay\arcs_overlay imi_' im_num_str];
        imwrite(arc_imi/255,[arcImPath '.png']);
        dlmwrite(arcs_write,bubble_info,'delimiter',',','newline','pc','Precision','%1.8e')
        my_bubble_stats=bubble_stats(vertex_total,bubble_info,arc_eqns,arc_errors);
        statsPath=[savePath '\stats\stats_bub_nxy_pace imi_' im_num_str];
        stats_write=[statsPath '_final.txt'];
        time_stamp=zeros(numel(my_bubble_stats(:,1)),1)+i1;
        dlmwrite(stats_write,[my_bubble_stats,time_stamp],'delimiter',',','newline','pc','Precision','%1.8e')
    end
end








