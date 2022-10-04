% +
% NAME: square_lattice_make
%
% PURPOSE:
%     This program will solve for h, the "hyperuniformity length" from the variance
%     of a point pattern.
%
% CATEGORY:
%     Hyperuniformity
%
% CALLING SEQUENCE:
%    find_h_length
%
% INPUTS: None at the moment as it is a stand alone program
%
% OPTIONAL INPUTS: (none)
%
% KEYWORD PARAMETERS: (none)
%
% OUTPUTS: Will write a text, that has h computed for different point patterns
%          read in

% SIDE EFFECTS: (none)
%
% MODIFICATION HISTORY:
%    written by: A. Chieco, UPenn, March  2018
%
%-
% clear all
% close all

num_bubs_vec=20*[1,2,3,5,10];
rescale=1;
xmin=0; xmax=rescale;
ymin=0; ymax=rescale;

for n_count=1:numel(num_bubs_vec)
    num_bubs=num_bubs_vec(n_count);
    for i0=1:1
        my_try=rescale*rand(num_bubs,2);        
        %we need to tag the vertice with each other. Thevx and vy will for
        %closed loops but do not give us information about which points lie
        %within the voronoi cell. v_pos and connections does tell us which
        %point is surrounded by which vertices but only has one vertex at
        %infinity,infiinity to describe all of the cell vertices outside
        %some bounded region. We identify the correct loops of vertices
        %around a point.
        [v_pos,connections]=voronoin(my_try);
        [vx,vy]=voronoi(my_try(:,1),my_try(:,2));
        [new_v_pos,new_cells]=voronoi_vertex_connect(v_pos,connections,vx,vy,my_try);
        
        hold off
        plot(my_try(:,1),my_try(:,2),'or');
        hold on
        plot(v_pos(:,1),v_pos(:,2),'ok')
        plot(new_v_pos(:,1),new_v_pos(:,2),'.m')
        vert_id_count=numel(new_v_pos(:,1));
        
        %we are going to construct a matrix with cell ID, poisson point and
        %list of ids of the vertices that surround the poisson point. Later
        %we will use these points and vertices to find centroids and areas
        %for the cells
        cells=zeros(num_bubs,20)-1;
        for i1=1:num_bubs
            cells(i1,1:3)=[i1,my_try(i1,:)];
            vert_ids=new_cells{i1};
            %Here we are going to adjust the ids of the vertices. We
            %may have to add new vertices to th cell of a bubble
            %because the current ones go over the edge. If that is the
            %cas we will add the new vertex to the vertex list, add an
            %id to the wertex list and switch that new ID into the ids
            %of the vertices for the cell
            for i2=1:numel(vert_ids)
                if i2<numel(vert_ids)
                    v1=new_v_pos(vert_ids(i2),:);
                    v1_place=i2;
                    v2=new_v_pos(vert_ids(i2+1),:);
                    v2_place=i2+1;
                else
                    v1=new_v_pos(vert_ids(i2),:);
                    v1_place=i2;
                    v2=new_v_pos(vert_ids(1),:);
                    v2_place=1;
                end
                if or(and(and(v1(1)>=xmin,v1(1)<=xmax),and(v1(2)>=xmin,v1(2)<=xmax)),...
                        and(and(v2(1)>=xmin,v2(1)<=xmax),and(v2(2)>=xmin,v2(2)<=xmax)))==1
                    slope=(v2(2)-v1(2))/(v2(1)-v1(1));
                    b1=v2(2)-slope*v2(1);
                    %we check the x values of the first vertex
                    if or(v1(1)<=0,v1(1)>rescale)==1
                        if v1(1)<=0
                            x_new=0;
                        else
                            x_new=rescale;
                        end
                        y_new=slope*x_new+b1;
                        %These check if we are off a diagnoal from the bounding box
                        if y_new>rescale
                            v1=[(rescale-b1)/slope,rescale];
                        elseif y_new<0
                            v1=[(rescale-b1)/slope,0];
                        else
                            %This is the standard
                            v1=[x_new,y_new];
                        end
                        vert_id_count=vert_id_count+1;
                        new_v_pos=[new_v_pos;[v1,vert_id_count]];
                        vert_ids(v1_place)=vert_id_count;
                    end
                    %we check the x values of the second vertex
                    if or(v2(1)<=0,v2(1)>rescale)==1
                        if v2(1)<=0
                            x_new=0;
                        else
                            x_new=rescale;
                        end
                        y_new=slope*x_new+b1;
                        %These check if we are off a diagnoal from the bounding box
                        if y_new>rescale
                            v2=[(rescale-b1)/slope,rescale];
                        elseif y_new<0
                            v2=[(rescale-b1)/slope,0];
                        else
                            %This is just the standard
                            v2=[x_new,y_new];
                        end
                        vert_id_count=vert_id_count+1;
                        new_v_pos=[new_v_pos;[v2,vert_id_count]];
                        vert_ids(v2_place)=vert_id_count;
                    end
                    %We check the y values if they are above the bounding box for
                    %the first vertex
                    if or(v1(2)<=0,v1(2)>rescale)==1
                        if v1(2)<=0
                            y_new=0;
                        else
                            y_new=rescale;
                        end
                        x_new=(y_new-b1)/slope;
                        %We don't need to check off the diagonal because this is
                        %covered in the x check. WE only do the standard value
                        v1=[x_new,y_new];
                        vert_id_count=vert_id_count+1;
                        new_v_pos=[new_v_pos;[v1,vert_id_count]];
                        vert_ids(v1_place)=vert_id_count;
                    end
                    %We check the y values if they are above the bounding box for
                    %the second vertex
                    if or(v2(2)<=0,v1(2)>rescale)==1
                        if v2(2)<=0
                            y_new=0;
                        else
                            y_new=rescale;
                        end
                        x_new=(y_new-b1)/slope;
                        %We don't need to check off the diagonal because this is
                        %covered in the x check. WE only do the standard value
                        v2=[x_new,y_new];
                        vert_id_count=vert_id_count+1;
                        new_v_pos=[new_v_pos;[v2,vert_id_count]];
                        vert_ids(v2_place)=vert_id_count;
                    end
                    dist=sqrt((v2(1)-v1(1))^2+(v2(2)-v1(2))^2);
                    if abs(v1(1)-v2(1))>abs(v1(2)-v2(2))
                        xs=(v1(1):(v2(1)-v1(1))/(1E4):v2(1));
                        ys=slope*xs+b1;
                    else
                        ys=(v1(2):(v2(2)-v1(2))/(1E4):v2(2));
                        xs=(ys-b1)/slope;
                    end
                    xy_mat=[xs',ys'];
                    %                 xy_keep1=xy_mat(and(xy_mat(:,1)>=0,xy_mat(:,1)<=rescale),:);
                    %                 xy_keep=xy_keep1(and(xy_keep1(:,2)>=0,xy_keep1(:,2)<=rescale),:);
                    %                 plot(xy_keep(:,1),xy_keep(:,2),'.b');
                    plot(xy_mat(:,1),xy_mat(:,2),'.b');
                    axis equal
                    %                 xlim([0 rescale])
                    %                 ylim([0 rescale])
                    bark=0;
                else
                    keyboard
                end
            end
            cells(i1,4:4+numel(vert_ids)-1)=vert_ids;
            i1
        end
        keyboard  
    end
end

%         vx=vx';
%         vy=vy';
%         
%         v1s=[vx(:,1),vy(:,1)];
%         v2s=[vx(:,2),vy(:,2)];
%         
% %         plot(vx,vy,'Linewidth',2,'Color','b')
%         hold on
%         plot(my_try(:,1),my_try(:,2),'or','Linewidth',2)
%         axis equal
%         xlim([0 rescale])
%         ylim([0 rescale])
%         voronoi(my_try(:,1),my_try(:,2));
%         
%         imi_mat=zeros(rescale,rescale);
%         
%         im_path='F:\Chieco\Hyperuniformity\HU foams\voronoi images\matlab write';
%         im_name=['\voronoi nbubs=' num2str(num_bubs) ' rescale=' num2str(rescale) ' imi_' num2str(i0) '.png'];        
%         imwrite(imi_mat/255,[im_path im_name],'png')
%         file_folder='F:\Chieco\Hyperuniformity\HU foams\voronoi images';
%         file_name=['\voronoi nbubs=' num2str(num_bubs) ' rescale=' num2str(rescale) ' imi_' num2str(i0) '.txt'];
%         dlmwrite([file_folder file_name],my_try,'Delimiter',',','Precision','%1.15e','newline','pc')










