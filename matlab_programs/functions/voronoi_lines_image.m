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

rescale=2848;
num_bubs_vec=3019*[1,2,3,5,10];

for n_count=1:numel(num_bubs_vec)
    num_bubs=num_bubs_vec(n_count);
    for i0=1:1       
        my_try=(rescale-1)*rand(num_bubs,2)+1;
        
        [vx,vy]=voronoi(my_try(:,1),my_try(:,2));
        vx=vx';
        vy=vy';
        
        v1s=[vx(:,1),vy(:,1)];
        v2s=[vx(:,2),vy(:,2)];
        
%         plot(vx,vy,'Linewidth',2,'Color','b')
        hold on
        plot(my_try(:,1),my_try(:,2),'or','Linewidth',2)
        axis equal
        xlim([0 rescale])
        ylim([0 rescale])
        voronoi(my_try(:,1),my_try(:,2));
        
        imi_mat=zeros(rescale,rescale);
        for i1=1:numel(vx(:,1))
            slope=(v2s(i1,2)-v1s(i1,2))/(v2s(i1,1)-v1s(i1,1));
            b1=v2s(i1,2)-slope*v2s(i1,1);
            %we make sure we are smapling enough points to get a good
            %voronoi diagram
            %we check the x values of the first vertex 
            if or(v1s(i1,1)<=0,v1s(i1,1)>rescale)==1
                if v1s(i1,1)<=0
                    x_new=0;
                else
                    x_new=rescale;
                end
                y_new=slope*x_new+b1;
                %These check if we are off a diagnoal from the bounding box
                if y_new>rescale
                    v1s(i1,:)=[(rescale-b1)/slope,rescale];
                elseif y_new<0
                    v1s(i1,:)=[(rescale-b1)/slope,0];
                else
                    %This is the standard
                    v1s(i1,:)=[x_new,y_new];
                end
            end
            %we check the x values of the second vertex 
            if or(v2s(i1,1)<=0,v2s(i1,1)>rescale)==1
                if v2s(i1,1)<=0
                    x_new=0;
                else
                    x_new=rescale;
                end
                y_new=slope*x_new+b1;
                %These check if we are off a diagnoal from the bounding box
                if y_new>rescale
                    v2s(i1,:)=[(rescale-b1)/slope,rescale];
                elseif y_new<0
                    v2s(i1,:)=[(rescale-b1)/slope,0];
                else
                    %This is just the standard
                    v2s(i1,:)=[x_new,y_new];
                end
            end
            %We check the y values if they are above the bounding box for
            %the first vertex
            if or(v1s(i1,2)<=0,v1s(i1,2)>rescale)==1
                if v1s(i1,2)<=0
                    y_new=0;
                else
                    y_new=rescale;
                end
                x_new=(y_new-b1)/slope;
                %We don't need to check off the diagonal because this is
                %covered in the x check. WE only do the standard value
                v1s(i1,:)=[x_new,y_new];
            end
            %We check the y values if they are above the bounding box for
            %the second vertex
            if or(v2s(i1,2)<=0,v1s(i1,2)>rescale)==1
                if v2s(i1,2)<=0
                    y_new=0;
                else
                    y_new=rescale;
                end
                x_new=(y_new-b1)/slope;
                %We don't need to check off the diagonal because this is
                %covered in the x check. WE only do the standard value
                v2s(i1,:)=[x_new,y_new];
            end
            dist=sqrt((v2s(i1,1)-v1s(i1,1))^2+(v2s(i1,2)-v1s(i1,2))^2);
            if abs(v1s(i1,1)-v2s(i1,1))>abs(v1s(i1,2)-v2s(i1,2))
                xs=(v1s(i1,1):(v2s(i1,1)-v1s(i1,1))/(1E4):v2s(i1,1));
                ys=slope*xs+b1;
            else
                ys=(v1s(i1,2):(v2s(i1,2)-v1s(i1,2))/(1E4):v2s(i1,2));
                xs=(ys-b1)/slope;
            end
            xy_mat=round([xs',ys']);
            xy_keep1=xy_mat(and(xy_mat(:,1)>=1,xy_mat(:,1)<=rescale),:);
            xy_keep=xy_keep1(and(xy_keep1(:,2)>=1,xy_keep1(:,2)<=rescale),:);
            for i2=1:numel(xy_keep(:,1))
                imi_mat(xy_keep(i2,2),xy_keep(i2,1))=200;
            end
        end        
        im_path='F:\Chieco\Hyperuniformity\HU foams\voronoi images\matlab write';
        im_name=['\voronoi nbubs=' num2str(num_bubs) ' rescale=' num2str(rescale) ' imi_' num2str(i0) '.png'];        
        imwrite(imi_mat/255,[im_path im_name],'png')
        file_folder='F:\Chieco\Hyperuniformity\HU foams\voronoi images';
        file_name=['\voronoi nbubs=' num2str(num_bubs) ' rescale=' num2str(rescale) ' imi_' num2str(i0) '.txt'];
        dlmwrite([file_folder file_name],my_try,'Delimiter',',','Precision','%1.15e','newline','pc')
    end
end












