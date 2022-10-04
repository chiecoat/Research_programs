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
var_read_path=[file_path '\real space averages\variance\foam data'];
% var_write_path=[file_path '\variance\list averages'];

hud_write_path=[file_path '\real space averages\hu length\foam data'];

data_type={'arc','vert'};
im_num_vec={'125','250','500'};
r_start=1; r_fin=1;

for imi_count=1:3
    im_num=im_num_vec{imi_count};
    for dat=1:2%numel(data_type)
        %first we read in each file
        fname=['E:\Chieco\Hyperuniformity\HU foams\patterns_to_analyze\'...
            data_type{dat} '_stats_mat imi_' im_num '.txt'];
        mat=dlmread(fname);
        if strcmp(data_type{dat},'arc')==1
            len_pad=sum(mat(:,3).^2)/sum(mat(:,3));
%             mat_keep_tot=mat(and(and(mat(:,1)>x_min+len_pad,mat(:,1)<x_max-len_pad),...
%                 and(mat(:,2)>y_min+len_pad,mat(:,2)<y_max-len_pad)),:);
%             L_sys_x=x_max-x_min-2*len_pad;
%             L_sys_y=y_max-y_min-2*len_pad;
            mat_keep_tot=mat(and(and(mat(:,1)>x_min,mat(:,1)<x_max),...
            and(mat(:,2)>y_min,mat(:,2)<y_max)),:);
            L_sys_x=x_max-x_min;
            L_sys_y=y_max-y_min;
            A_sys=L_sys_x*L_sys_y;
            L_sys=min([L_sys_x,L_sys_y]);
            sys_specs=[x_min,y_min;x_max,y_max];
        else
            mat_keep_tot=mat(and(and(mat(:,1)>x_min,mat(:,1)<x_max),...
                and(mat(:,2)>y_min,mat(:,2)<y_max)),:);
            L_sys_x=x_max-x_min;
            L_sys_y=y_max-y_min;
            A_sys=L_sys_x*L_sys_y;
            L_sys=min([L_sys_x,L_sys_y]);
            sys_specs=[x_min,y_min;x_max,y_max];
        end
        for runs=1:1%numel(mat(1,3:end))+1
            if runs==1
                areas=mat_keep_tot(:,3);
                type='points_weighted';
                phi=sum(mat_keep_tot(:,3))/A_sys;
            elseif runs==numel(mat(1,3:end))+1
                areas_mono=mat_keep_tot(:,3)*0+1;
                type='points_mono';
            else
                continue
                %areas=mat(:,end);
                %type='centers_coord';
            end
            variance_path=[var_read_path '\variance_circle_windows_avgellnorm_interiorOnly ' data_type{dat} '_' type ' imi_' im_num '_avg.txt'];
            variance_data=dlmread(variance_path);
            figure(1)
            loglog(2*variance_data(:,1),variance_data(:,6),'Linewidth',2)
            hold on;
            %plot(2*variance_data(:,1),particle_area./(pi*(variance_data(:,1)*length_norm).^2),'--r','LineWidth',3);
            %if and(runs==1,w==numel(weight_vec))
            %   xs=logspace(-4,4);
            %   plot(xs,4./(pi*xs.^2),'--r','LineWidth',3);
            %end
            length_norm=variance_data(1,2);
            particle_area=sum((mat_keep_tot(:,3)).^2)/sum(mat_keep_tot(:,3));%variance_data(1,3);
            %some data is held didfferenty, this is just so I can change them to whatever format
            %I need and not have to change any variables
            sqrt_avga=sqrt(particle_area(1));
            if strcmp('points_weighted',type)==1
                if strcmp(data_type{dat},'bub')==1
                    norm_for_vals=sqrt_avga;
                else
                    norm_for_vals=sum(mat_keep_tot(:,3).^2)/sum(mat_keep_tot(:,3));
                end
                bub_length_scale=sqrt_avga;
                norm_var=sum(areas)/A_sys;
                use_for_h=particle_area(1);
            else
                norm_var=1;
                use_for_h=numel(areas)/A_sys;
                norm_for_vals=sqrt(A_sys/numel(areas));
            end
            %use this for the voronoi packings
            if runs==1
                variance_in=[variance_data(:,1)*length_norm,variance_data(:,6),zeros(numel(variance_data(:,1)),1)];
            else
                variance_in=[variance_data(:,1)*length_norm,variance_data(:,6),zeros(numel(variance_data(:,1)),1)];
            end
            err_calc=zeros(numel(variance_in(:,1)),4);
            %Now we need to calculate the error for the variance in phi
            %for continuum model we have S=[1-(1-q)^M ]/q   where  q=Vwin/Vsys
            for err=1:numel(err_calc(:,1))
                qs=(variance_in(err,1)^2)/A_sys;
                S_num=(1-(1-qs)^variance_data(err,5))/qs;
                err_small=(variance_data(err,4)/variance_in(err,1)^6)/(S_num*norm_var);
                err_big=(2*variance_in(err,2)^2)/(S_num-1);
                err_calc(err,1)=sqrt(err_small+err_big);
                err_calc(err,2:4)=[qs,S_num,variance_in(err,1)/sqrt_avga];
            end
            variance_in(:,3)=err_calc(:,1);
            if runs==1
                [hs,error]=calculate_h_circle(variance_in,use_for_h);
            else
                [hs,error]=calculate_h_circle_points(variance_in,use_for_h);
            end
            figure(2)
%             errorbar(2*hs(:,1)/norm_for_vals,hs(:,2)/norm_for_vals,real(error)/norm_for_vals,'LineWidth',2)
            errorbar(2*hs(:,1)/norm_for_vals,hs(:,2)/norm_for_vals,real(error)/norm_for_vals,'LineWidth',2)
            set(gca,'YScale','log','XScale', 'log')
            %         h_look=[mean(hs(2*hs(:,1)/norm_for_vals>6,2)),std(hs(2*hs(:,1)/norm_for_vals>6,2))]/norm_for_vals;
            hold on
            
            hud_file_write=[hud_write_path '\hu length circle_windows ' data_type{dat} '_' type ' imi_' im_num '_avg.txt'];
            dlmwrite(hud_file_write,[2*hs(:,1)/norm_for_vals,hs(:,2)/norm_for_vals,real(error)/norm_for_vals],'Newline','pc');
        end
    end
end
loglog(2*hs(:,1)/norm_for_vals,hs(:,1)/norm_for_vals,'--r','LineWidth',3)
xs_small=logspace(-4,1,100);
ys_small=xs_small/2.*(1-sqrt(phi*norm_for_vals*(pi/4)*(xs_small.^2)));
plot(xs_small,ys_small,'Linewidth',2)

