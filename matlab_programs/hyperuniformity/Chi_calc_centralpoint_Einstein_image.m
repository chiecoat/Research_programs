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
% clear all;
% close all;

sys_size_vec=[2^12];

for my_sz=1:numel(sys_size_vec)
    %There are three different types of packings
    file_path='E:\Chieco\Hyperuniformity\HU defects\area defects\lattice triangular data\N=1E6';
    num_in=11547624;
    num_in_txt='11547624';
    r_start=1; r_fin=5;    
    %These are the area fractions we want to investigate
    p_f_vec=90;
    f_path=[file_path '\fourier space data'];    
    for t1=1:numel(p_f_vec)
        p_frac=p_f_vec(t1);
        txt_phi=num2str(p_f_vec(t1));
        %first we read in each file
        for read=r_start:r_fin             
            deets='delta_0.5b f=0 ad_mag=0';
            fname=[file_path '\particle positions\' ...
                'Einstein_triangle_defect N' num_in_txt ' phi_' txt_phi '_loop_' num2str(read) ' ' deets '.txt'];
            mat=dlmread(fname);
            %We have and exact expression for the average area or a particle when we
            %are in the continuum limit. We compute the area here, knowing the ratio
            %between the two species of particles is 1.4
            x_min=0; x_max=mat(2,1); y_min=0; y_max=mat(2,2);
            sys_specs=[x_min,y_min;x_max,y_max];
            L_sys_x=x_max-x_min;
            L_sys_y=y_max-y_min;
            A_sys=L_sys_x*L_sys_y;
            mat=mat(3:end,:);
            lists=1;%numel(mat(1,:))/3;
            N_tot=numel(mat(:,1));
            N_big=N_tot;
            if read==r_start
                phi=p_frac/10^(ceil(log10(p_frac+1)));
                norm=phi;
                %Now we compute diameter of a small particle, sig_s
                sig_s=sqrt(4*phi*A_sys/(N_tot*pi));
                a_small=pi*(sig_s/2)^2;
                %here is the the phi weighted average area of a particle, <a>
                phi_small=N_tot*a_small/A_sys;                               
            end
            avg_a=sum(mat(:,3).^2)/sum(mat(:,3)); 
            r_vec=sqrt(mat(:,3)/pi);
            edges_cut=20*sqrt(avg_a);
            %We need two vectors for each particle population.
            sz_scale_x=ceil((L_sys_x+edges_cut)*sys_size_vec(my_sz));
            sz_scale_y=ceil((L_sys_y+edges_cut)*sys_size_vec(my_sz));
            rescale=min([sz_scale_y,sz_scale_x]);
            sz_area=sz_scale_x*sz_scale_y;
            edges_rescale=ceil(edges_cut*rescale/2);
            for my_lists=0:lists-1
                mat_in=rescale*[double(mat(:,3*my_lists+1)-x_min),double(mat(:,3*my_lists+2)-y_min),r_vec];
                %we rescale our average a
                phi_small_scale=sum(pi*min(mat_in(:,3)).^2)/sz_area;
                avg_a_scale=pi*sum(mat_in(:,3).^4)/sum(mat_in(:,3).^2);
                mat_bi=[ceil(mat_in(:,1)),ceil(mat_in(:,2)),pi*mat_in(:,3).^2];
                mat_bi(mat_bi(:,1)==0,1)=zeros(numel(mat_bi(mat_bi(:,1)==0,1)),1)+1;
                mat_bi(mat_bi(:,2)==0,2)=zeros(numel(mat_bi(mat_bi(:,2)==0,2)),1)+1;  
                clear mat
                clear mat_in
%                 my_inds=sub2ind([sz_scale_y,sz_scale_x],mat_bi(:,2),mat_bi(:,1));
%                 rep=histc(my_inds,(1:sz_scale_y*sz_scale_x));
%                 imi_in=zeros(rescale,rescale);
%                 imi_in(1:end)=rep*avg_a_scale;
                imi_in=sparse(mat_bi(:,2),mat_bi(:,1),mat_bi(:,3),sz_scale_y,sz_scale_x);                
                imi_in=imi_in(edges_rescale:rescale-edges_rescale,edges_rescale:rescale-edges_rescale);
                size(imi_in)
                imi_norm=1/(sqrt(sum(sum((imi_in).^2))));
                spec_add=fft_radial_large_edit(imi_in,imi_norm,0.1,sqrt(avg_a_scale));
                clear imi_in
                clear mat_bi
                spec_tot=[spec_add(spec_add(:,3)>0,1),spec_add(spec_add(:,3)>0,3)];
                loglog(spec_tot(2:end,1),spec_tot(2:end,2),'Linewidth',2)%'Color', [255, 215, 0]./255)
                hold on
                bark=0;
                bidi_name=['centralpoint Chi_q_try_Einstein_' deets ' phi_' txt_phi '_' num2str(read) ' N_sys_' num_in_txt ' L_sys_' num2str(sz_scale_x) '.txt'];
                dlmwrite([f_path '\' bidi_name],spec_add,'Precision',12,'newline','pc');
            end            
        end
    end
end