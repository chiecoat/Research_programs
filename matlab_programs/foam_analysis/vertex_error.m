% +
% NAME: vertex_error
%
% PURPOSE:%     
%
% CATEGORY:
%     Foam coarsening
%
% CALLING SEQUENCE:
%    vert_loc_error=vert_loc_error(vert_mat,foam_imi)
%
% INPUTS: 
%    vert_,at: this is a Nx8 matrix where the columns are
%    id,x,y,theta, location 'energy', bubble IDs (all will have 3)
%    [id_vert,x_vert,y_vert,theta,energy,cell_id_1,cell_id_2,cell_id_3] 
%
%    foam_imi: this is a NxM image of the foam. The matrix x,y locations
%    of the vertices correspond with the column,row indices of the image
%
%    box_size: we need an LxL box to test for the local energy landscape.
%    This box has to be at least 3 and must be odd
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
%    written by: A. Chieco UPenn, September 2017
%-
function [vert_loc_error]=vertex_error(vert_mat,vert_stack,mask,foam_imi,box_size,theta_check)


%We know the locations of the vertices but want to test the space around
%them for error in that location. We will assign energy values to x and y
%positions within a small square surrounding the "known" positions. This
%will wieght each surrounding positins
n_verts=numel(vert_mat(:,1));
sz_vert=numel(vert_stack(:,1,1));
sz_half=floor(sz_vert/2);
num_x=numel(foam_imi(1,:));
num_y=numel(foam_imi(:,1));
x_check_vec=(-floor(box_size/2):floor(box_size/2));
y_check_vec=(-floor(box_size/2):floor(box_size/2));
vert_loc_error=[vert_mat,zeros(numel(vert_mat(:,1)),3)];


for i1=1:n_verts
    %this is the condition for a boundary vertex
    if vert_mat(i1,4)<0
        continue
    end
    %This is the x,y,theta positions of the vertex
    loc=vert_mat(i1,2:4);
    loc_E=vert_mat(i1,5);
    energy_mat=zeros(box_size,box_size);
    x_mat=energy_mat;
    y_mat=energy_mat;
    %we got across x positions
    for x_i=1:box_size
        %we also check along the y positions
        for y_i=1:box_size
            if and(x_i==0,y_i==0)==1
                x_mat(y_i,x_i)=loc(1);
                y_mat(y_i,x_i)=loc(2);
                energy_mat(y_i,x_i)=loc_E;
            else
                row=loc(2)+y_check_vec(y_i);
                col=loc(1)+x_check_vec(x_i);
                if or(or(row-sz_half<1,row+sz_half>num_y),or(col-sz_half<1,col+sz_half>num_x))==1
                    x_mat(y_i,x_i)=col;
                    y_mat(y_i,x_i)=row;
                    energy_mat(y_i,x_i)=1E7;                    
                else
                    x_mat(y_i,x_i)=col;
                    y_mat(y_i,x_i)=row;
                    foam_mat=foam_imi(row-sz_half:row+sz_half,col-sz_half:col+sz_half).*mask;
                    energy_mat(y_i,x_i)=sum(sum((foam_mat-vert_stack(:,:,loc(3))).^2))/(sz_vert^2);
                end
            end
        end
    end
    weight_mat=1./((energy_mat).^2);
    x_bar=sum(sum(x_mat.*weight_mat))/(sum(sum(weight_mat)));
    y_bar=sum(sum(y_mat.*weight_mat))/(sum(sum(weight_mat)));
    x_err=2*abs(loc(1)-x_bar);
    y_err=2*abs(loc(2)-y_bar);
    vert_loc_error(i1,end-2:end-1)=[x_err,y_err];
    %Once we know the error in the location of the matrix we want to know
    %the error in its angle. We will now take the location as known and
    %rotate the vertex through some angles. We will fit a parabola to the
    %energy versus theta landscape and the minimum will be the actual
    %theta. The error in the fit will be the error in theta.
    theta_mat=zeros(2*theta_check+1,2);
    theta_vec_keep=1;
    for t_i=-theta_check:theta_check;
        theta=loc(3)+t_i;
        if theta>120
            theta_mark=theta-120;
       %Theta is either 0 or negative so we have to add 120 to
        elseif theta<=0
            theta_mark=120+theta;
        elseif and(theta>=1,theta<=120)==1
            theta_mark=theta;
        end
        foam_mat=foam_imi(loc(2)-sz_half:loc(2)+sz_half,loc(1)-sz_half:loc(1)+sz_half).*mask;
        energy_loc=sum(sum((foam_mat-vert_stack(:,:,theta_mark)).^2))/(sz_vert^2);
        theta_mat(theta_vec_keep,:)=[theta,energy_loc];     
        theta_vec_keep=theta_vec_keep+1;
    end
%     plot(theta_mat(:,1),theta_mat(:,2));    
    theta_coeffs=polyfit(theta_mat(:,1),theta_mat(:,2),2);
    theta_new=-theta_coeffs(2)/(2*theta_coeffs(1));
    theta_err=2*abs(loc(3)-theta_new);
    %This means that we could not fit a parabola well to the location.
    %However we know that this is the minimum theta location, so the error
    %will just be the total degrees we scanned.
    if theta_err>theta_check
       theta_err=theta_check;
    end
    vert_loc_error(i1,end)=[theta_err];    
end

if numel(vert_loc_error((vert_loc_error(:,4)<0),:))>0
    vert_loc_error((vert_loc_error(:,4)<0),end-2)=max(vert_loc_error(:,end-2));
    vert_loc_error((vert_loc_error(:,4)<0),end-1)=max(vert_loc_error(:,end-1));
    vert_loc_error((vert_loc_error(:,4)<0),end)=max(vert_loc_error(:,end));
end
