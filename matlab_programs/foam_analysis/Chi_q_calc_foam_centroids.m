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
% clear all 
% close all

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
%     rgb = label2rgb(Ld,'jet',[.5 .5 .5]);
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
    %This is for the foam centroids
    for i1=1:numel(hu_mat(:,1))
        imi_centroids(hu_mat(i1,2),hu_mat(i1,1))=hu_mat(i1,3);
    end
    
    avg_a=sum(areas.^2)/(A_sys);
    avg_a_cubed=sum(areas.^4)/(A_sys);
    N_tot=numel(hu_mat(:,1));
    
    a_nums=unique(imi_centroids);

    imi_norm=1/sqrt(sum(a_nums.^2));
    imi_in=sparse(imi_centroids);
    clear imi_centroids
    spec_add=fft_radial_large(imi_in,imi_norm);
    clear imi_in
    spec_tot=[spec_add(:,1)*sqrt(avg_a)/sqrt(A_sys),spec_add(:,2)];
%     loglog(spec_tot(2:end,1),spec_tot(2:end,2))%'Color', [255, 215, 0]./255)
%     hold on
    bidi_name=['fftshift_binary Structure_f poly imi_' num2str(num_foam) '.txt'];
    dlmwrite([f_path '\' bidi_name],spec_tot,'Precision',12,'newline','pc');
    bark=0;
end

%     close Figure 1