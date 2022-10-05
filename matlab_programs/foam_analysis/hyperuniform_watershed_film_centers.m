% +
% NAME: hyperuniform_watershed_film_centers
%
% PURPOSE:%   
%     This program will convert a skeletonized image of films into a point
%     pattern, using the midpoint of the films as the point in the pattern.
%     This 
%
% CATEGORY:
%     Hyperuniformity
%
% CALLING SEQUENCE:
%     hyperuniform_watershed_film_centers
%
% INPUTS: 
%    None, but we will need to know what image files to read in.
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
% clear all 
% close all

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
   %we find the x and y of foam frame, plot the x and y in a figure
    thresh = find(Ld==0);
    xy_mat=[-1,-1];
    for i1=1:numel(Ld(1,:));
        ys=find(Ld(:,i1)==0);
        xy_new=zeros(numel(ys),2);
        xy_new(:,1)=i1;
        xy_new(:,2)=ys;
        xy_mat=[xy_mat;xy_new];
    end
    xy_mat=xy_mat(2:end,:);
%     plot(xy_mat(:,1),xy_mat(:,2),'.r');
%     hold on

   %here we find the location of the foam vertices
    xy_tag=[xy_mat,zeros(numel(xy_mat(:,1)),3)];
    xy_tag(:,5)=-1;
    for i2=1:numel(xy_mat(:,1))
        % establish a small box named Ld_small at each x and y
        if or(xy_mat(i2,2)-1<1,xy_mat(i2,2)+1>numel(Ld(:,1)))==1
            continue
        elseif or(xy_mat(i2,1)-2<1,xy_mat(i2,1)+1>numel(Ld(1,:)))==1
            continue
        end
        Ld_small=Ld((xy_mat(i2,2)-1):(xy_mat(i2,2)+1),(xy_mat(i2,1)-1):(xy_mat(i2,1)+1));
       % find all the different values in the small box
        unique_vals=unique(Ld_small);
        % find the non-zero value
        bub_id=find(unique_vals~=0);
        if numel(bub_id)>3
            keyboard
        end
        % add the non-zero values into third through n columns of i2 row
        xy_tag(i2,3:3+numel(bub_id)-1)=unique_vals(bub_id);
    end
    % keep the first and second columns of xy_tag if the value of fifth column is greater than
    % zero
    xy_keep=xy_tag(xy_tag(:,5)>0,:);
%     plot(xy_keep(:,1),xy_keep(:,2),'og')
    
    L_sys_x=numel(Ld(1,:));
    L_sys_y=numel(Ld(:,1));
    A_sys=L_sys_x*L_sys_y;   
    %We make a binary image for the film skeleton
    imi_films=zeros(L_sys_y,L_sys_x);
    imi_films(find(Ld==0))=1;
    imi_films_poisson=make_poisson_films(imi_films,xy_keep);    
    
    [film_centers,film_lengths]=find_film_centers(imi_films,xy_keep);
%     plot(film_centers(:,3),film_centers(:,4),'ob')
%     axis equal
    
    hu_mat=[film_centers(:,3:4),film_lengths];
    hu_mat(hu_mat(:,1)<1,1)=1;
    hu_mat(hu_mat(:,2)<1,2)=1;
    hu_mat(hu_mat(:,1)>L_sys_x,1)=L_sys_x;
    hu_mat(hu_mat(:,2)>L_sys_y,2)=L_sys_y;   
    
    %Now we either search for HU in centroids or in the binary foam image
    imi_film_centers=zeros(L_sys_y,L_sys_x);
    %This is for the foam centroids
    for i1=1:numel(hu_mat(:,1))
        imi_film_centers(hu_mat(i1,2),hu_mat(i1,1))=hu_mat(i1,3);
    end
    
    avg_ell=sum(film_lengths.^2)/sum(film_lengths);
    avg_ell_cubed=sum(film_lengths.^4)/sum(film_lengths);

    for i2=1:1
        r_start=1; r_fin=1;
        if i2==1
            imi=imi_film_centers;
            phi=sum(sum(imi_film_centers))/A_sys;
            norm=phi
            L_norm=sqrt(avg_ell)
            L_norm^2
            continue
%             [L_norm,sqrt(2*phi/(3*numel(xy_keep(:,1))))*L_sys_x]
            tag='film_centers_binomial';
            imwrite(imi,[imiPath '\imi_' num2str(imi_nums(num_foam)) ' film_centers.png'],'png')
        end  
        if i2==2
            imi=imi_film_centers;
            imi(imi_film_centers>0)=1;
            phi=sum(sum(imi))/A_sys;
            norm=phi;
            L_norm=1;            
            tag='film_points_binomial';
            imwrite(imi,[imiPath '\imi_' num2str(imi_nums(num_foam)) ' film_centers.png'],'png')
        end 
        %we find how many windows and hat values here
        c_rads1=logspace(log10(1E-4*sqrt(avg_ell)),log10(sqrt(avg_ell)),51);
        c_rads2=logspace(log10(sqrt(avg_ell)),log10(0.8*L_sys),151);
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
                mat_write_stack(i3,:,i4)=[L_in/L_norm,avg_ell,avg_ell_cubed,n_in,var_matlab_loop/norm];
            end
            %Now we compute the 'corrected two pass variance' from Numerical
            %recipes. It is also called 'compensated variant' on wikipedia.
            %This is what Matlab uses to calculate the variance
            var_matlab=var(phi_vec);
            vec_write(i3,:)=[L_in/L_norm,avg_ell,avg_ell_cubed,n_in,var_matlab/norm];
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
% 
% 
% 
% 
% 
% 
