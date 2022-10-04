% +
% NAME: Chi_q_calc_sparse_image
%
% PURPOSE:
%     This program will find the spectral density for bidisperse mixture of
%     disks. To do so it will call the "partial structure" code. Patrial
%     structure can be found in the "Function" folder in MatLab Programs.
%
% CATEGORY:
%     Hyperuniformity
%
% CALLING SEQUENCE:
%    FFT_my_image
%
% INPUTS: 
%    There are no inputs but we need to read in particle positoins for
%    whatever configuration we are trying to analyze.
%
% OPTIONAL INPUTS: (none)
%
% KEYWORD PARAMETERS: (none)
%
% OUTPUTS: Chi_q
%    The function will write data files that are the spectral density 
%    versus wave number for our input configuration.     
%
% 
% SIDE EFFECTS: (none)
%
% MODIFICATION HISTORY:
%    written by: Remi Dreyfus, UPenn
%     edited by: A. Chieco, UPenn 2017
%-
clear all 
close all

savePath='F:\Chieco\Hyperuniformity\HU foams';
imiPath=[savePath '\images'];
imi_nums=[130,337,530];

f_path=[savePath '\fourier space data']

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
%     plot(xy_keep(:,1),xy_keep(:,2),'ob')
    
    L_sys_x=numel(Ld(1,:));
    L_sys_y=numel(Ld(:,1));
    A_sys=L_sys_x*L_sys_y;   
    %We make a binary image for the film skeleton
    imi_films=zeros(L_sys_y,L_sys_x);
    imi_films(find(Ld==0))=1;
    imi_films_poisson=make_poisson_films(imi_films,xy_keep);    
    
    film_lengths=find_film_lengths(imi_films,xy_keep);
    
    avg_ell=sum(film_lengths.^2)/sum(film_lengths);

    imi_norm=1/sqrt(sum(film_lengths.^2));
    imi_in=sparse(imi_films);
    clear imi_centroids
    spec_add=fft_radial_large(imi_in,imi_norm);
    clear imi_in
    spec_tot=[spec_add(:,1)*sqrt(avg_ell)/sqrt(A_sys),spec_add(:,2)];
    loglog(spec_tot(2:end,1),spec_tot(2:end,2))%'Color', [255, 215, 0]./255)
    hold on
    bidi_name=['fftshift_binary Structure_f poly_films imi_' num2str(imi_nums(num_foam)) '.txt'];
    dlmwrite([f_path '\' bidi_name],spec_tot,'Precision',12,'newline','pc');
    bark=0;
end