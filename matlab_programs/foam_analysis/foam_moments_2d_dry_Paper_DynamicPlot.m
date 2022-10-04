% +
% NAME: foam_track_2d
%
% PURPOSE:%
%
% CATEGORY:
%     Foam coarsening
%
% CALLING SEQUENCE:
%    foam_track_2d
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
%    written by: A. Chieco, UPenn, December 2017
% % This is the folder used
% clear all
% close all

date_vec={'2021-09-23','2021-11-23','2021-11-29','2021-12-01'};

num_start=0;
loop_start_vec=[1,24,45,40];
loop_end_vec=[69,124,124,200];

color_mat=[0,0,204;...
    204,102,0;...
    0,104,0;...
    204,0,204]/255;

p_color_mat=[102,102,102;...
             204,0,0;...
             255,140,0;...
             51,255,0;...
             0,153,153;...
             0,0,204;...
             204,51,255]/255;

% date_vec={'2021-11-30'};
% num_start=0;
% loop_start_vec=[1];
% loop_end_vec=[60];
symbol_vec={'o','+','x'};
m_size=[8,12,12];


for dates=1:3%numel(date_vec)
    folder_titles=date_vec{dates};
    %First we read in all of the files with the bubbles information
    imiPath=['E:\foam coarsening\dry foam\' folder_titles '\binary edges'];
    filePath=['E:\Chieco\Dry Foam Moments\' folder_titles '\bubble data\'];
    
    num_start=0;
    loop_start=loop_start_vec(dates);
    loop_end=loop_end_vec(dates);
    
    %The conversion from pixel to mm is 1mm=22.005pix  and we convert to real
    %units in the when reading in the files.
    mm_conv=38.9;
    
    %first we plot time sequence data fr <<n>>_p and <A^p>/<A>^p
    p_max=21;
    powers=(1:p_max);
    % n_times=numel(times);
    n_powers=numel(powers);
    
    N_bub_mat=zeros(1,4);
    N_bub_mat_plot=zeros(1,4);
    
    num_step=1;
    list_length=numel((loop_start:num_step:loop_end));
    
    A_moment_raw_data=zeros(list_length,n_powers);
    avg_avg_n_raw_data=zeros(list_length,n_powers);
    
    k_o_calc_vec=[0.0252,0.0228,0.0223];
    
    for i1=loop_start:num_step:loop_end
        file=[filePath '\bubbles imi_' num2str(num_start+i1) '_moments.txt'];
        %The bubble information we read in is
        %n_sides,x_bub_cen,y_bub_cen,perimeter,del_perimeter,area,del_area,
        %circularity,del_circularity,elongation,del_elongation, "time". The del
        %indicates uncertainty in that measurment The "time" is just how many
        %frames we are from the initial frame
        bubble_info=dlmread(file);
        %In order for track to work we need to have the coordinates of particle
        %centers first
        bub_stack_all=[bubble_info(:,1:3),bubble_info(:,4)./(mm_conv^2),bubble_info(:,5:6),bubble_info(:,7)/mm_conv,zeros(numel(bubble_info(:,1)),1)+i1];
        bub_stack_all=bub_stack_all(and(bub_stack_all(:,4)>0,bub_stack_all(:,4)<1E4),:);
        bub_stack=bub_stack_all;
        
        min_area=0;
        
        sides_min=3;
        sides_max=300;
        bub_stack=bub_stack(and(bub_stack(:,3)>=sides_min,bub_stack(:,3)<=sides_max),:);
        bub_stack=bub_stack(bub_stack(:,4)>min_area,:);
        
        %we need to identify the particle centers
        x_cens=bub_stack(:,1);
        y_cens=bub_stack(:,2);
        areas=bub_stack(:,4);
        circ=bub_stack(:,5);
        ecc=bub_stack(:,6);
        perim=bub_stack(:,7);
        bub_stack=bub_stack(bub_stack(:,end)>0,:);
        
        %plot(bub_stack(:,1),bub_stack(:,2),'ok','Linewidth',2)
        
        n_sides_vec=(2:20)';%unique(bub_stack(:,3));
        areas=bub_stack(:,4);
        
        
        area_tot=sum(areas);
        avg_a=mean(areas);
        p_start=-5;
        p_end=21;
        powers=(p_start:1:p_end);
        for i2=1:numel(powers)
            power=powers(i2);
            avg_a_power=mean(areas.^power);
            A_moment_raw_data(i1-loop_start+1,i2)=avg_a_power;
            %now we need to calculate <<n>>_p
            n_data_vec=zeros(1,numel(n_sides_vec));
            if i2==1
                F_of_n_mat=zeros(numel(n_sides_vec),2);
            end
            for i3=1:numel(n_sides_vec)
                n_sides=n_sides_vec(i3);
                n_sides_data=bub_stack(bub_stack(:,3)==n_sides,:);
                if numel(n_sides_data)==0
                    continue
                end
                n_data_vec(i3)=(n_sides*sum(n_sides_data(:,4).^power))/(sum(areas.^power));
                if i2==1
                    F_of_n_mat(i3,:)=[n_sides,sum(n_sides_data(:,4))/area_tot];
                end
            end
            avg_avg_n_raw_data(i1-loop_start+1,i2)=sum(n_data_vec);
            %         	plot(F_of_n_mat(:,1),F_of_n_mat(:,2))
        end
    end
    
    %we have to shift the times so <A(t=0)>=0
    times=transpose(5*((1:numel(A_moment_raw_data(:,1)))-1));
    avg_a_vs_time_plot=[times,A_moment_raw_data(:,7)];
    area_fit_lin=@(coeffs,times) coeffs(1)*(times-coeffs(2));
    coeffs_in=[0.025,100];
    [coeffs_lin,Rs,Jac,CovB,MSE]=nlinfit(avg_a_vs_time_plot(:,1),avg_a_vs_time_plot(:,2),area_fit_lin,coeffs_in);
    
    times_shift=times-coeffs_lin(2);
    
    
    
% %     p_color_mat=hsv(numel(p_vals));
%     
%     x_weights_for_vals=zeros(1,numel(col_vals));
%     %this is so we can plot a line for K_o behind the data
%     i_for_ko=1;
%     p_val=p_vals(i_for_ko);
%     col_val=i_for_ko+6;
%     y_vals=A_moment_raw_data(:,col_val+1)./A_moment_raw_data(:,col_val);
%     avg_avg_n_of_p=mean(avg_avg_n_raw_data(:,col_val));
%     x_weight_for_vals=(1+1/p_val)*(avg_avg_n_of_p-6);
%     hold on
%     k_o
    p_vals=[1,2,3,5,7,10,15];
    col_vals=p_vals+6;
    for fit_loop=1:2
        data_tot=[0,0];
        for i4=1:numel(col_vals)
            column=col_vals(i4);
            p_val=p_vals(i4);
            y_vals=A_moment_raw_data(:,column+1)./A_moment_raw_data(:,column);
            avg_avg_n_of_p=mean(avg_avg_n_raw_data(:,column));
            x_weights_for_vals(:,i4)=(1+1/p_val)*(avg_avg_n_of_p-6);
            x_weight_for_vals=(1+1/p_val)*(avg_avg_n_of_p-6);
            if fit_loop==1
                data_tot=[data_tot;x_weight_for_vals*times_shift,y_vals];                
            else
                plot(x_weight_for_vals*times_shift,y_vals,symbol_vec{dates},'Linewidth',2,'MarkerSize',12,'Color',p_color_mat(i4,:))
                hold on
            end
        end
        data_tot=data_tot(2:end,:);
        if and(fit_loop==1,dates==2)==1
            all_data_A_B=data_tot;
        elseif and(fit_loop==1,dates==3)==1
            all_data_A_B=[all_data_A_B;data_tot];
        end
        t_vals=(0:50:6E3);
        if fit_loop==1
            [data,indices]=sort(data_tot(:,1));
            k_o=polyfit(data_tot(indices,1),data_tot(indices,2),1);
            plot(t_vals,t_vals*k_o(1)+k_o(2),'--k','Linewidth',2)
            hold on
            plot(t_vals,t_vals*k_o_calc_vec(dates),'--g','Linewidth',2)
            k_o
        end
        if and(fit_loop==1,dates==3)==1
            [data,indices]=sort(all_data_A_B(:,1));
            k_o_combo=polyfit(all_data_A_B(indices,1),all_data_A_B(indices,2),1);
            plot(t_vals,t_vals*k_o_combo(1)+k_o_combo(2),'--k','Linewidth',2)
            hold on
            k_o_combo
        end
    end
    
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.03],'LineWidth',3,'FontSize',20)
end

