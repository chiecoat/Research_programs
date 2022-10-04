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

file_path='E:\Chieco\Hyperuniformity\HU foams';
data_type={'vert','arc','bub'};
r_start=1; r_fin=1;

%the normalization length for all of the data is just the sqrt(<a>) of the
%bubble.
bub_num='500';
foam_norm=dlmread(['E:\Chieco\Hyperuniformity\HU foams'...
                  '\patterns_to_analyze\bub_stats_mat square_imi_' bub_num '.txt']);
area_norm=sum(foam_norm(:,3).^2)/sum(foam_norm(:,3));
bub_length_norm=sqrt(area_norm);

weight_vec=[1];

for dat=3:3%numel(data_type)
    f_path=[file_path '\fourier space data'];
    %first we read in each file
    fname=[file_path  '\patterns_to_analyze\'...
        data_type{dat} '_stats_mat square_imi_' bub_num '.txt'];
    mat=dlmread(fname);
    if dat==1
        x_mid=(min(mat(:,1))+max(mat(:,1)))/2;
        y_mid=(min(mat(:,2))+max(mat(:,2)))/2;
        L_sys_x=max(mat(:,1))-min(mat(:,1));
        L_sys_y=max(mat(:,2))-min(mat(:,2));
        L_sys=min([L_sys_x,L_sys_y]);
        if L_sys==L_sys_x
            x_min=min(mat(:,1)); x_max=max(mat(:,1));
            y_min=y_mid-L_sys/2; y_max=y_mid+L_sys/2-1;
        else
            x_min=x_mid-L_sys/2; x_max=x_mid+L_sys/2-1;
            y_min=min(mat(:,2)); y_max=max(mat(:,2));
        end
        sz_scale=ceil(L_sys);
        A_sys=L_sys^2;
    end
    for w=1:numel(weight_vec)
        width=weight_vec(w);
        for runs=1:1%numel(mat(1,3:end))+1
            if runs==1
                mat=mat(and(and(mat(:,1)>x_min,mat(:,1)<x_max),...
                       and(mat(:,2)>y_min,mat(:,2)<y_max)),:);
                areas=mat(:,3);
                type='points_weighted';
            elseif runs==numel(mat(1,3:end))+1
                areas=mat(:,3)*0+1;
                type='points_mono';
            else
                continue
                %             areas=mat(:,end);
                %             type='centers_coord';
            end
            a_vec=width*areas;
            if runs==1
                if strcmp(data_type{dat},'bub')==1
                    length_norm=sqrt(sum(areas.^2)/sum(areas));
                else
                    length_norm=sum(a_vec.^2)/sum(a_vec);
                end
            else
                length_norm=sqrt(A_sys/numel(areas));
            end
%             length_norm=bub_length_norm;
            mat_in_tot=[mat(:,1:2),a_vec];
            mat_in_x=mat_in_tot(and(mat_in_tot(:,1)>x_min,mat_in_tot(:,1)<x_max),:);
            mat_in=mat_in_x(and(mat_in_x(:,2)>y_min,mat_in_x(:,2)<y_max),:);
            %We fill in the matrix
            mat_bi=[ceil(mat_in(:,1)-x_min)+1,ceil(mat_in(:,2)-y_min)+1,mat_in(:,3)];
            mat_bi(mat_bi(:,1)>sz_scale,1)=sz_scale;
            mat_bi(mat_bi(:,2)>sz_scale,2)=sz_scale;
            imi_norm=1/sqrt(sum(mat_bi(:,3).^2));
            imi_in=sparse(mat_bi(:,2),mat_bi(:,1),mat_bi(:,3),sz_scale,sz_scale);
            spec_add=fft_radial_large_edit(imi_in,imi_norm,10,0.1,length_norm);
            clear imi_in
            clear mat_bi
            spec_tot=[spec_add(spec_add(:,3)>0,1)/(2*pi),spec_add(spec_add(:,3)>0,3)];
            min_check=[length_norm,min(spec_tot(spec_tot(:,1)>0,1))]
            loglog(spec_tot(:,1),spec_tot(:,2),'Linewidth',2)%'Color', [255, 215, 0]./255)
            hold on
            bark=0;
            keep=isnan(spec_add(:,3));
            spec_write=[spec_add(keep==0,1:2)/(2*pi),spec_add(keep==0,3:4)];
            bidi_name=['\Chi_q_avgellnorm ' data_type{dat} '_' type ' width_eq_' num2str(width) ' square_imi_' bub_num ' L_sys_' num2str(L_sys) '.txt'];
            dlmwrite([f_path bidi_name],spec_write,'Precision',12,'newline','pc');            
        end
    end
end

