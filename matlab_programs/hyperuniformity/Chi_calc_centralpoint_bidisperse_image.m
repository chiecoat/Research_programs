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
tic
sys_size_vec=[2^11,2^12];

for my_sz=1:numel(sys_size_vec)
    %There are three different types of packings
    for titles=2:2
        if titles==1
            file_path='E:\Chieco\Hyperuniformity\HU Jammed Packings from Carl Goodrich';
            jam_type='\quenchFromTinf';
            tol_vec={''};
            num_in=2048;
            num_in_txt='2048';
            r_start=1; r_fin=1;
            %We have and exact expression for the average area or a particle when we
            %are in the continuum limit. We compute the area here, knowing the ratio
            %between the two species of particles is 1.4
            x_min=0; x_max=1; y_min=0; y_max=1;
            L_sys=x_max-x_min; sig_con=1.4;
            A_sys=L_sys^2;
            quench_rate={''};
            %These are the area fractions we want to investigate
            p_f_vec=[8400];%[(1000:500:7500),(8000:100:9000),(8420:20:8480)];
        end
        if titles==2
            file_path='E:\Chieco\Hyperuniformity\HU Jammed Packings from Ning Xu';
            jam_type='\quenchFromTinf';
            tol_vec={'\tol1.0E-12'};
            num_in=1E5;
            num_in_txt='1E5';
            r_start=1; r_fin=5;
            %We have and exact expression for the average area or a particle when we
            %are in the continuum limit. We compute the area here, knowing the ratio
            %between the two species of particles is 1.4
            x_min=-0.5; x_max=0.5; y_min=-0.5; y_max=0.5;
            L_sys=x_max-x_min; sig_con=1.4;
            A_sys=L_sys^2;
            quench_rate={''};
            %These are the area fractions we want to investigate
            p_f_vec=[8412];
        end
        if titles==3
            file_path='E:\Chieco\Hyperuniformity\HU Jammed Packings from Ning Xu';
            jam_type='\quenchVaryRate';
            tol_vec={'\tol1.0E-12\N=1E6'};
            num_in=1E6;
            num_in_txt='1E6';
            r_start=1; r_fin=1;
            %We have and exact expression for the average area or a particle when we
            %are in the continuum limit. We compute the area here, knowing the ratio
            %between the two species of particles is 1.4
            x_min=-0.5; x_max=0.5; y_min=-0.5; y_max=0.5;
            L_sys=x_max-x_min; sig_con=1.4;
            A_sys=L_sys^2;
            quench_rate={'step_0'};
            %These are the area fractions we want to investigate
            p_f_vec=[85,86];
        end
        f_path=[file_path jam_type tol_vec{1} '\fourier space data'];
        sys_specs=[x_min,y_min;x_max,y_max];
        for t1=1:numel(p_f_vec)
            tol=tol_vec{1};
            p_frac=p_f_vec(t1);
            txt_phi=num2str(p_f_vec(t1));
            for p_f=1:numel(quench_rate)
                %first we read in each file
                for read=r_start:r_fin
                    if titles==1
                        fname=[file_path jam_type tol  '\particle positions\' ...
                            '\0p' txt_phi 'xyrMats_200.txt'];
                        mat=dlmread(fname);
                        lists=1;%numel(mat(1,:))/3;
                    end
                    
                    if titles==2
                        fname=[file_path jam_type tol  '\particle positions\config2d\'...
                              'config2d_100000_0.8412_2.0_1.4_' num2str(read) '.txt'];
%                             '\phi.' txt_phi '_' num2str(read)];
                        mat=dlmread(fname);
                        mat=mat(2:end,:);
                        lists=1;
                    end
                    if titles==3
                        fname=[file_path jam_type tol  '\particle positions\CSV\' quench_rate{p_f} ' CSV'...
                            '\confire_0.' txt_phi '_' num2str(read)];
%                         mat=dlmread(fname,',',[1 0 num_in 1]);
                        mat=dlmread(fname);
                        lists=1;
                    end
                    N_tot=numel(mat(:,1));
                    N_big=N_tot/2; N_small=N_tot/2;
                    if read==r_start
                        phi=p_frac/10^(ceil(log10(p_frac+1)));
                        norm=phi;
                        %Now we compute diameter of a small particle, sig_s
                        sig_s=sqrt(8*phi*A_sys/(N_tot*pi*(1+sig_con^2)));
                        a_small=pi*(sig_s/2)^2;
                        a_big=pi*(sig_con*sig_s/2)^2;
                        %here is the the phi weighted average area of a particle, <a>
                        phi_big=N_big*a_big/A_sys; phi_small=N_small*a_small/A_sys;
                        avg_a=(phi_big*a_big+phi_small*a_small)/phi;
                        r_vec=[zeros(N_tot/2,1)+sig_s/2;zeros(N_tot/2,1)+sig_con*sig_s/2];
                    end
                    %We need two vectors for each particle population.
                    sz_scale=sys_size_vec(my_sz);
                    for my_lists=0:lists-1
                        mat_in=sz_scale*[double(mat(:,3*my_lists+1)-x_min),double(mat(:,3*my_lists+2)-y_min),r_vec];
                        %we rescale our average a
                        phi_big_scale=N_big*pi*max(mat_in(:,3)).^2/sz_scale^2;
                        phi_small_scale=N_small*pi*min(mat_in(:,3)).^2/sz_scale^2;
                        avg_a_scale=(phi_big_scale*pi*max(mat_in(:,3)).^2+phi_small_scale*pi*min(mat_in(:,3)).^2)/phi;
                        mat_bi=[ceil(mat_in(:,1)),ceil(mat_in(:,2)),pi*mat_in(:,3).^2];
                        mat_bi(mat_bi(:,1)==0,1)=zeros(1,numel(mat_bi(mat_bi(:,1)==0,1)))+1;
                        mat_bi(mat_bi(:,2)==0,2)=zeros(1,numel(mat_bi(mat_bi(:,2)==0,2)))+1;
                        imi_norm=1/sqrt(N_small*(pi*min(mat_in(:,3)).^2)^2+N_big*(pi*max(mat_in(:,3)).^2)^2);
%                         if titles>1
%                             clear mat
%                             clear r_vec
%                         end
%                         clear mat_in                      
                        imi_in=sparse(mat_bi(:,2),mat_bi(:,1),mat_bi(:,3),sz_scale,sz_scale);
                        dq=0.1;
                        spec_add=fft_radial_large_edit(imi_in,imi_norm,dq,sqrt(avg_a_scale));
                        clear imi_in
                        spec_tot=[spec_add(spec_add(:,3)>0,1),spec_add(spec_add(:,3)>0,3)];
                        loglog(spec_tot(2:end,1)/(2*pi),spec_tot(2:end,2),'Linewidth',2)%'Color', [255, 215, 0]./255)
                        hold on
                        bark=0;
%                         bidi_name=['centralpoint Chi_q_log_' num2str(dq) ' ' quench_rate{p_f} 'phi_' txt_phi '_' num2str(my_lists+1) ' N_sys_' num_in_txt ' L_sys_' num2str(sz_scale) '.txt'];
%                         dlmwrite([f_path '\' bidi_name],spec_add,'Precision',12,'newline','pc');
                    end
                end
            end
        end
    end
end
toc