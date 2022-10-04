% +
% NAME: watershed_bubbles
%
% PURPOSE:%     
%
% CATEGORY:
%     Foam coarsening
%
% CALLING SEQUENCE:
%    vertex_mat_tot=zoom_vertex_connect(vert_mat)
%
% INPUTS: 
%    vert_mat: this is a Nx8 matrix where the columns are
%    id,x,y,theta, vertex type, location 'energy'
%    [id_vert,x_vert,y_vert,theta,energy] 
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
%    written by: A. Chieco UPenn, October 2018
%-
function [vertex_mat_tot]=zoom_vertex_connect(vertex_total,sz_vert)

num_verts=numel(vertex_total(:,1));
mask_rad=sz_vert/2;

vertex_mat_tot=[vertex_total,zeros(num_verts,3)-1];
% a=imread('F:\foam coarsening\20181018 zoom_test\DSC images\DSC_3128.jpg');
% color=[1,0,0;0,1,0;0,1,1];

for i1=1:num_verts
    %%he "primary leg" which was tracked in the vertex finding code will be
    %used to make a fiduciary line that seeds our vertex connectiong code.
    %The angle that is tracke dis in realtion to treating the center of the
    %structuring element as (0,0)
    theta_0=-vertex_total(i1,4)*pi/180;
    %***We use a negative theta for the y values because of a quirk in how the
    %image is represented. the (0,0) for the image is the top left but the
    %angle was calculated in reference to a (0,0) in the bottom right of
    %the image. This is annoying but shouldn't affect the analysis
    xy_edge=[vertex_total(i1,2)+mask_rad*cos(theta_0),vertex_total(i1,3)+mask_rad*sin(theta_0)];
    %We need to find the "inverse" locations of the vertex legs. This will
    %define a region where we look for neighboring vertices
    inv_angles=(-vertex_total(i1,4)+([60,-60,-180]))*pi/180;
    inv_angles(inv_angles<0)=inv_angles(inv_angles<0)+2*pi;
    inv_angles=sort(inv_angles);
    %We translate all the points to polar coordinates with the vertex 
    %whose neighbors we are trying to find as the (0,0) location for the 
    %Cartesian points
    [theta,r_s]=cart2pol(vertex_total(:,2)-vertex_total(i1,2),vertex_total(:,3)-vertex_total(i1,3));
    theta(theta<0)=theta(theta<0)+2*pi;
    vert_info_total=[vertex_total(:,1:3),theta,r_s];
    [dist_mins,dist_ind]=sort(vert_info_total(:,5));
    vert_info_check=vert_info_total(dist_ind,:);   
    theta=theta(dist_ind);
%     imshow(a);
%     hold on
    for i2=1:3 
        if i2==1
            vert_check=vert_info_check(and(theta>=inv_angles(1),theta<=inv_angles(2)),:);
        elseif i2==2
            vert_check=vert_info_check(and(theta>=inv_angles(2),theta<=inv_angles(3)),:);
        else
            vert_check=vert_info_check(or(theta>=inv_angles(3),theta<=inv_angles(1)),:);
        end
        vert_keep=vert_check(vert_check(:,1)~=i1,:);        
        if numel(vert_keep)==0
            continue
        end
        true_neighbor=0;
        for i3=1:numel(vert_keep(:,1))
            a_sq=sum((vertex_total(i1,2:3)-xy_edge(1:2)).^2);
            b_sq=sum((vertex_total(i1,2:3)-vert_keep(i3,2:3)).^2);
            c_sq=sum((xy_edge(1:2)-vert_keep(i3,2:3)).^2);
            denom=-2*sqrt(a_sq)*sqrt(b_sq);
            law_cos_theta=acos((c_sq-a_sq-b_sq)/denom)*180/pi;
            if or(and(law_cos_theta>90,law_cos_theta<150),law_cos_theta<30)==1
                v_neighbor_id=vert_keep(i3);
                true_neighbor=1;
                break
            end
        end
        %This means we have identified a real neighbor and not jsut
        %vertices that exist in the large swath we are searching in for
        %neighbors
        if true_neighbor==1
            v_neighbor_loc=find(vertex_mat_tot(i1,end-2:end)<0)-1;
            vertex_mat_tot(i1,end-2+v_neighbor_loc(1))=v_neighbor_id;            
        end                
    end
end

%Now that we have connected all the vertices we want to check if the
%connections are true. Wrong connections are made because vertex A finds
%vertex B but vertex B does not find vertex A. 
for i4=1:num_verts
    v_id=vertex_mat_tot(i4,1);
    neighbor_ids=vertex_mat_tot(i4,end-2:end);
    n_ids=neighbor_ids(neighbor_ids>0);
    count=0;
    for i5=1:numel(n_ids)
        neighbor_check=vertex_mat_tot(n_ids(i5),end-2:end);
        %if both vertices have each other as neighbors then we do not have
        %to do anything. However if this is zero that means only one vertex
        %has identified the other as a neighbor
        if numel(neighbor_check(neighbor_check==v_id))==0
            n_ids(i5)=-1;
            count=count+1;
        end
%         else
%            plot([vertex_mat_tot(i4,2),vertex_mat_tot(n_ids(i5),2)],...
%                 [vertex_mat_tot(i4,3),vertex_mat_tot(n_ids(i5),3)],...
%                 '-o','Linewidth',4)
%             hold on
%         end            
    end
    if count>0
        neighbor_ids_replace=n_ids(n_ids>0);
        vertex_mat_tot(i4,end-2:end)=[-1,-1,-1];
        vertex_mat_tot(i4,end-2:end-2+(numel(neighbor_ids_replace)-1))=neighbor_ids_replace;
    end
    
end











end