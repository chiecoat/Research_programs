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
%    written by: A. Chieco & Y. Shen, UPenn, September 2017
%-
clear all 
close all

savePath='F:\Chieco\Hyperuniformity\HU foams';
imiPath=[savePath '\images'];
imi_nums=[130,337,530];

for num_foam=1:numel(imi_nums)    
    vertexPath=[imiPath '\imi_' num2str(imi_nums(num_foam)) ' square_crop_skeleton.png'];
    vertex_imi=imread(vertexPath);
    L_sys=numel(vertex_imi(:,1));
    
    vertex_imi(vertex_imi<150)=0;
    vertex_imi(vertex_imi>150)=255;
    
    extra1=3;
    watershed_imi=zeros(L_sys+2*extra1,L_sys+2*extra1);
    %we are going to pad watershed imi and then draw a box around the foam
    watershed_imi(extra1:end-extra1-1,extra1:end-extra1-1)=vertex_imi;
    watershed_imi(1:extra1-1,1:end)=255;
    watershed_imi(end-extra1:end,1:end)=255;
    watershed_imi(1:end,1:extra1-1)=255;
    watershed_imi(1:end,end-extra1:end)=255;
    
    L_watershed=numel(watershed_imi(:,1));
    extra2=2*extra1;
    watershed_full=zeros(L_watershed+2*extra2,L_watershed+2*extra2);
    %we are going to pad watershed imi and then draw a box around the foam
    watershed_full(extra2:end-extra2-1,extra2:end-extra2-1)=watershed_imi;
    
    Ld_full=watershed(watershed_full);
    Ld=Ld_full(extra2+extra1:end-extra2-extra1-1,extra2+extra1:end-extra2-extra1-1);
    rgb = label2rgb(Ld,'jet',[.5 .5 .5]);
%     imshow(rgb);
    % we find the x and y of foam frame, plot the x and y in a figure
    thresh = find(Ld==0);
    xy_mat=[-1,-1];
    
    stats=regionprops(Ld,'Centroid','Area');
    areas=cat(1, stats.Area);
    centroids=cat(1, stats.Centroid);
    
    L_sys_x=numel(Ld(1,:));
    L_sys_y=numel(Ld(:,1));
    A_sys=L_sys_x*L_sys_y;
    
    hu_mat=[centroids,areas];
    hu_mat=[round(hu_mat(hu_mat(:,3)>0,1)),round(hu_mat(hu_mat(:,3)>0,2)),hu_mat(hu_mat(:,3)>0,3)];
    hu_mat(hu_mat(:,1)<1,1)=1;
    hu_mat(hu_mat(:,2)<1,2)=1;
    hu_mat(hu_mat(:,1)>L_sys_x,1)=L_sys_x;
    hu_mat(hu_mat(:,2)>L_sys_y,2)=L_sys_y;   
    
    %Now we either search for HU in centroids or in the binary foam image
    imi_centroids=zeros(L_sys_y,L_sys_x);
    imi_films=imi_centroids;
    %This is for the foam centroids
    for i1=1:numel(hu_mat(:,1))
        imi_centroids(hu_mat(i1,2),hu_mat(i1,1))=hu_mat(i1,3);
    end
    %This is for foam skeleton
    imi_films(find(Ld==0))=1;
    numel(centroids(:,1))
    continue
    
    for i2=1:1
        r_start=1; r_fin=1;
        if i2==1
            avg_a=sum(areas.^2)/sum(areas);
            avg_a_cubed=sum(areas.^4)/sum(areas);
            N_tot=numel(hu_mat(:,1));
            a_vec=hu_mat(:,3);
            phi=1-sum(sum(imi_films))/A_sys;
            norm=phi;
            L_norm=sqrt(avg_a);
            continue
            imi=imi_centroids;            
            tag='centroids';
        else
            avg_a=1;
            avg_a_cubed=1;
            N_tot=numel(hu_mat(:,1));
            a_vec=zeros(numel(hu_mat(:,3)),1)+1;  
            imi=imi_centroids;
            imi(imi_centroids>0)=1;
            phi=sum(sum(imi))/A_sys;
            norm=phi;
            L_norm=1;                      
            tag='centroids_points';
        end
        %we find how many windows and hat values here
        c_rads1=logspace(log10(1E-4*sqrt(avg_a)),log10(sqrt(avg_a)),50);
        c_rads2=logspace(log10(sqrt(avg_a)),log10(0.8*L_sys),151);
        c_rads=ceil([c_rads1,c_rads2(2:end)]);
        c_rads=unique(c_rads);
        N_windows_vec=zeros(1,numel(c_rads))+1E5;
        %Now we test for hyperuniformity at a given box size L and keep all of
        %the raw data in a matrix
        i1_start=1; i1_fin=numel(c_rads);
        vec_write=zeros(i1_fin,5);
        mat_write_stack=zeros(i1_fin,5,r_fin);
        for i3=1:i1_fin
            n_in=N_windows_vec(i3);
            L_in=c_rads(i3);
            n_phi_vec=r_fin*n_in;
            phi_vec=zeros(n_phi_vec,1);
            phi_sq_vec=phi_vec;
            for i4=r_start:r_fin
                phi_list=square_window_pixel(imi,n_in,L_in);
                mat_start=(i4-1)*n_in+1; mat_fin=i4*n_in;
                phi_vec(mat_start:mat_fin,1)=phi_list(:,1);
                phi_sq_vec(mat_start:mat_fin,1)=phi_list.^2;
                var_matlab_loop=var(phi_list);
                mat_write_stack(i3,:,i4)=[L_in/L_norm,avg_a,avg_a_cubed,n_in,var_matlab_loop/norm];
            end
            %Now we compute the 'corrected two pass variance' from Numerical
            %recipes. It is also called 'compensated variant' on wikipedia.
            %This is what Matlab uses to calculate the variance
            var_matlab=var(phi_vec);
            vec_write(i3,:)=[L_in/L_norm,avg_a,avg_a_cubed,n_in,var_matlab/norm];
            text_path1=[savePath '\real space averages\variance'];
            filename=['\variance ' tag ' imi_' num2str(imi_nums(num_foam)) ' square_crop_skeleton.txt'];
            dlmwrite(strcat(text_path1,filename),vec_write(1:i3,:),'Delimiter',',','Precision','%1.15e','newline','pc');
        end
%         text_path2=strcat(file_path,jam_type,tol,'\loop data\variance');
%         for i3=r_start:r_fin
%             var_mat_write=mat_write_stack(:,:,i3);
%             filename=['\variance exact_circle_window_overlap ' txt_phi ' loop_',num2str(i3),'.txt'];
%             dlmwrite(strcat(text_path2,filename),var_mat_write,'Delimiter',',','Precision','%1.15e','newline','pc');
%         end
    end
 
    
end
