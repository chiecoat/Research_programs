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
%    written by: A. Chieco, UPenn, June 2022
% % This is the folder used
% clear all
% close all

date_vec={'2022-09-26 nitrogen too mono','2022-09-05 nitrogen scaling','2022-09-13 nitrogen too poly'};
num_start=0;

%this is the time is took to set up the foam
t_offset=[90,0,120];

%In date 1 the data sets are separated by 13 hours and 3 minutes because of
%a program malfunction. I am trying to stitch together the data sets
% nums_date_1=[(1:6:637),(647:6:1241)];
% times_date_1=[5*((1:6:637)-1),13*60+3+5*((647:6:1241)-1)]+t_offset(1);
nums_date_1=(48:6:1938);
times_date_1=5*nums_date_1+t_offset(1);

nums_date_2=[(72:6:1188),(1440:6:1741),(1194:6:1434)];
times_date_2=[5*((72:6:1188)-1),21*60+5*((1440:6:1741)-1440)+5*1187,46*60+15+5*((1194:6:1434)-1194)+5*1187]+t_offset(2);
% times_date_2=[5*((72:6:1188)-1),,45*60+45+5*((1194:6:1434)-1194+1188)]+t_offset(2);

% nums_date_2=(72:6:1740);
% times_date_2=5*(nums_date_2-1)+t_offset(2);


% nums_date_3=[(517:6:637),(647:6:1241)];
% times_date_3=[5*((517:6:637)-1),13*60+3+5*((647:6:1241)-1)]+t_offset(3);

nums_date_3=[(48:6:1848)];
times_date_3=[5*(nums_date_3-1)]+t_offset(3);


color_mat=[0,0,204;...
    104,104,104;...
    204,0,0;
    204,102,0]/255;

color_mat_Rc=[0,102,204;...
    0,0,0;...
    204,0,204;
    204,0,0]/255;

symbol_vec={'o','+','x','o'};
m_size=[8,12,12];

for dates=1:3%numel(date_vec)
    folder_titles=date_vec{dates};
    %First we read in all of the files with the bubbles information
    filePath=['E:\Chieco\non-scaling state\' folder_titles '\bubble data'];
    
    %The conversion from pixel to mm is different depending on the date of the runs
    %becasue the experimental set up changed. We keep track here.
    if dates==1
        mm_conv=30;
        nums_date=nums_date_1;
        times_date=times_date_1;
    elseif dates==2
        mm_conv=30;
        nums_date=nums_date_2;
        times_date=times_date_2;
    elseif dates==3
        mm_conv=30;
        nums_date=nums_date_3;
        times_date=times_date_3;
    else
        mm_conv=30;
        nums_date=nums_date_4;
        times_date=times_date_4;
    end
    
    num_start=0;  
    
    %first we plot time sequence data fr <<n>>_p and <A^p>/<A>^p
    p_max=20;
    powers=(1:p_max);
    % n_times=numel(times);
    n_powers=numel(powers);
    
    N_bub_mat=zeros(1,6);
    N_bub_mat_plot=zeros(1,4);
    
    list_length=numel(nums_date);
    
    A_moment_raw_data=zeros(list_length,n_powers);
    avg_avg_n_raw_data=zeros(list_length,n_powers);
    
    color=hsv(list_length);
    color_count=0;
    
    for i1=1:list_length
        color_count=color_count+1;
        num=nums_date(i1);
        file=[filePath '\bubbles imi_' num2str(num) '_xyNA_CirEccPer_data_full.txt'];


            

        bubble_info=dlmread(file);
        %In order for track to work we need to have the coordinates of particle
        %centers first
        bub_stack_all=[bubble_info(:,1:3),bubble_info(:,4)./(mm_conv^2),bubble_info(:,5:6),bubble_info(:,7)/mm_conv,zeros(numel(bubble_info(:,1)),1)+i1];
%         bub_stack_all=bub_stack_all(and(bub_stack_all(:,4)>0,bub_stack_all(:,4)<1E4),:);
        bub_stack=bub_stack_all;
        
        min_area=0/(mm_conv^2);
        
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
        
        area_tot=sum(areas);
        avg_a=mean(areas);
        %now we need to calculate F(n)
        F_of_n_mat=zeros(numel(n_sides_vec),2);
        for i3=1:numel(n_sides_vec)
            n_sides=n_sides_vec(i3);
            n_sides_data=bub_stack(bub_stack(:,3)==n_sides,:);
            if numel(n_sides_data)==0
                continue
            end
            F_of_n_mat(i3,:)=[n_sides,sum(n_sides_data(:,4))/area_tot];
        end
        
        N_bub_mat=[N_bub_mat;[bub_stack(:,1)*0+i1,bub_stack(:,3:4),sort(bub_stack(:,4))/mean(bub_stack(:,4)),bub_stack(:,1:2)]];     
        bub_ids=(0:numel(areas)-1)';
        plot_vals=[sort(bub_stack(:,4)),bub_ids/(numel(areas)-1)];
        plot_keep=plot_vals;
                
%         %adding individual curves for the area distribution
%         figure(1)
%         semilogy(plot_keep(:,1)/(mean(plot_keep(:,1))),1-plot_keep(:,2),'LineWidth',1,'Color',color(color_count,:))
%         set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',3,'fontweight','bold','fontsize',16)
%         xlim([0 8])
%         ylim([2E-4 1])
%         hold on
%         %adding individual curves for the side number distribution
%         figure(2)
%         loop_bin_edges=(0.5:1:max(bub_stack(:,3))+0.5);
%         if i1==list_length
%             histogram(bub_stack(:,3),loop_bin_edges,'Normalization','pdf','EdgeColor','k','FaceColor','none','LineWidth',2,'DisplayStyle','stairs')    
%         else
%             histogram(bub_stack(:,3),loop_bin_edges,'Normalization','pdf','EdgeColor',color(color_count,:),'FaceColor','none','LineWidth',1,'DisplayStyle','stairs')
%         end
%         hold on
        N_bub_mat_plot=[N_bub_mat_plot;[bub_stack(:,1)*0+i1,bub_stack(:,3:4),bub_stack(:,4)/mean(bub_stack(:,4))]];     
    end
    
     if dates==1
        t_scaling=2500;%2000;
        t_approach_to_scaling=2300;
     elseif dates==2
        t_scaling=300;
        t_approach_to_scaling=0;
     elseif dates==3
         t_scaling=300;
         t_approach_to_scaling=0;
    end
    
    N_bub_mat=N_bub_mat(2:end,:);
    N_bub_mat_plot=N_bub_mat_plot(2:end,:);
    
    areas_all=N_bub_mat(:,3);
    bub_ids=(0:numel(areas_all)-1)';
    plot_vals=[sort(areas_all),bub_ids/(numel(areas_all)-1)];
    plot_keep=plot_vals;
    
    areas_all_plot=N_bub_mat_plot(:,3);
    bub_ids_plot=(0:numel(areas_all_plot)-1)';
    plot_vals_plot=[sort(N_bub_mat_plot(:,4)),bub_ids_plot/(numel(areas_all_plot)-1)];
    plot_keep_plot=plot_vals_plot;
    
    
    %these are the usual measures of the scaling state. The average area
    %should increase linearly with time and the other three values should
    %be constant 
    figure(3)
    times=unique(N_bub_mat(:,1));
    bin_bars=(2:20);
    %this is to find the average area
    avg_a_vs_time=zeros(numel(times),3);
    avg_a_vs_time(:,1)=times_date;
    %this is to find the average first moment
    avg_asq_over_avg_a_sq_vs_time=zeros(numel(times),2);
    avg_asq_over_avg_a_sq_vs_time(:,1)=times_date;
    %this is to find the average sides number
    avg_n_vs_time=zeros(numel(times),2);
    avg_n_vs_time(:,1)=times_date;
    %this is to find the average area weighted side number
    avg_avg_n_vs_time=zeros(numel(times),2);
    avg_avg_n_vs_time(:,1)=times_date;
    for t1=1:numel(times)
        F_of_n_vec=bin_bars*0;
        areas_time=N_bub_mat(N_bub_mat(:,1)==times(t1),:);
        avg_a_vs_time(t1,2)=mean(areas_time(:,3));
        avg_a_vs_time(t1,3)=mean(areas_time(areas_time(:,2)==6,3));
        avg_asq_over_avg_a_sq_vs_time(t1,2)=mean(areas_time(:,3).^2)/(mean(areas_time(:,3))^2);
        avg_n_vs_time(t1,2)=mean(areas_time(:,2));
        area_tot=sum(areas_time(:,3));
        for F_calc=1:numel(bin_bars)-1
            n_sides=bin_bars(F_calc);
            area_vec=areas_time(areas_time(:,2)==n_sides,3);
            if numel(area_vec)==0
                continue
            else
                F_of_n_vec(F_calc)=sum(area_vec)/area_tot;
            end
        end
        avg_avg_n_vs_time(t1,2)=sum(bin_bars.*F_of_n_vec);
    end
    
    n_bub_tot=numel(N_bub_mat(:,1))
    
    avg_a_vs_time_plot=avg_a_vs_time(avg_a_vs_time(:,2)>0,:);
    
    figure(3)
    %these are the average radii vs time
    rads_mean=sqrt(avg_a_vs_time(:,2)/pi);
    %these are the average critical radii vs time
    rad_c_zero_mean=sqrt(avg_a_vs_time(:,3)/pi);
    b_over_a_fit=@(coeffs,rs) coeffs(1)*(rs-coeffs(2))+coeffs(3);
    coeffs_in=[0.99,rads_mean(1),rad_c_zero_mean(1)];
    [coeffs_bovera,Rs,Jac,CovB,MSE]=nlinfit(rads_mean,rad_c_zero_mean,b_over_a_fit,coeffs_in);
    
    plot(rads_mean,rad_c_zero_mean,'o','MarkerSize',8,'Linewidth',2,'Color',color_mat(dates,:))
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',2,'fontweight','bold','fontsize',16)
    hold on
    xs_plot=(0:max(rads_mean)/100:2*max(rads_mean));
    plot(xs_plot,coeffs_bovera(1)*(xs_plot-coeffs_bovera(2))+coeffs_bovera(3),'--k')
    b_over_a=coeffs_bovera(1)
    r_zero_vec=[coeffs_bovera(2),coeffs_bovera(2)/rads_mean(1)]
    r_czero_vec=[coeffs_bovera(3),coeffs_bovera(3)/rad_c_zero_mean(1)]
    r_czero_o_r_zero=r_czero_vec(1)/r_zero_vec(1)
    
    figure(4)
    plot(avg_a_vs_time_plot(:,1),rad_c_zero_mean,'o','Linewidth',2,'MarkerSize',8,'Color',color_mat(dates,:))
    hold on
    plot(avg_a_vs_time_plot(:,1),rads_mean,'x','Linewidth',2,'MarkerSize',8,'Color',color_mat_Rc(dates,:))

    
    
    area_fit_lin=@(coeffs,times) coeffs(1)*(times-coeffs(2));
    coeffs_in=[0.025,100];
    [coeffs_lin,Rs,Jac,CovB,MSE]=nlinfit(avg_a_vs_time_plot(:,1),avg_a_vs_time_plot(:,2),area_fit_lin,coeffs_in);
    
    figure(5)
    plot(avg_a_vs_time_plot(:,1),avg_a_vs_time_plot(:,2),symbol_vec{dates},'Linewidth',2,'MarkerSize',8,'Color',color_mat(dates,:))
    hold on
    plot(avg_a_vs_time_plot(:,1),avg_a_vs_time_plot(:,3),symbol_vec{dates},'Linewidth',2,'MarkerSize',8,'Color',color_mat_Rc(dates,:))
%     loglog(avg_a_vs_time_plot(:,1),sqrt(avg_a_vs_time_plot(:,2)/pi)-r_zero_vec(1),symbol_vec{dates},'Linewidth',2,'MarkerSize',8,'Color',color_mat(dates,:))
%     loglog(avg_a_vs_time_plot(:,1),sqrt(avg_a_vs_time_plot(:,3)/pi)-r_czero_vec(1),symbol_vec{dates},'Linewidth',2,'MarkerSize',8,'Color',color_mat(dates+1,:))
    lin_xs=(0:5E3);
    late_t_fit_data=[avg_a_vs_time_plot(:,1),avg_a_vs_time_plot(:,2)];
    late_t_fit_data=late_t_fit_data(late_t_fit_data(:,1)>t_scaling,:);
    late_t_lin_coeffs=polyfit(late_t_fit_data(:,1),late_t_fit_data(:,2),1);
    
    plot(lin_xs,lin_xs*late_t_lin_coeffs(1)+late_t_lin_coeffs(2),'-.k','Linewidth',3)
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',2,'fontweight','bold','fontsize',16)
    da_dt_slope_show=[coeffs_lin(1),late_t_lin_coeffs]
    
    figure(6)
    plot(avg_avg_n_vs_time(:,1),avg_avg_n_vs_time(:,end),symbol_vec{dates},'Linewidth',2,'MarkerSize',8,'Color',color_mat(dates,:))
    hold on
    mean_avg_avg_n_vs_time=mean(avg_avg_n_vs_time(avg_avg_n_vs_time(:,1)>t_scaling,end));
    plot([0,1E4],[mean_avg_avg_n_vs_time,mean_avg_avg_n_vs_time],':','Linewidth',2,'Color',color_mat(dates,:));
        
    plot(avg_asq_over_avg_a_sq_vs_time(:,1),avg_asq_over_avg_a_sq_vs_time(:,2),symbol_vec{dates},'Linewidth',2,'MarkerSize',8,'Color',color_mat(dates,:))
    mean_avg_asq_over_avg_a_sq_vs_time=mean(avg_asq_over_avg_a_sq_vs_time(avg_asq_over_avg_a_sq_vs_time(:,1)>t_scaling,end));
    plot([0,1E4],[mean_avg_asq_over_avg_a_sq_vs_time,mean_avg_asq_over_avg_a_sq_vs_time],'-.','Linewidth',2,'Color',color_mat(dates,:));
    
    plot(avg_n_vs_time(:,1),avg_n_vs_time(:,end),symbol_vec{dates},'Linewidth',2,'MarkerSize',8,'Color',color_mat(dates,:))
    mean_avg_n_vs_time=mean(avg_n_vs_time(avg_n_vs_time(:,1)>t_scaling,end));
    plot([0,1E4],[mean_avg_n_vs_time,mean_avg_n_vs_time],'--','Linewidth',2,'Color',color_mat(dates,:));
    
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',2,'fontweight','bold','fontsize',16)
    
    k_0=coeffs_lin(1)*mean(avg_asq_over_avg_a_sq_vs_time(:,end))/(2*(mean(avg_avg_n_vs_time(:,end))-6));
    area_moment=mean(avg_asq_over_avg_a_sq_vs_time(:,end));
    area_moment_err=std(avg_asq_over_avg_a_sq_vs_time(:,end));
    avg_avg_n=mean(avg_avg_n_vs_time(:,end));
    avg_avg_n_err=std(avg_avg_n_vs_time(:,end));
    
    areas_plot=[plot_keep_plot(:,1)/mean(plot_keep_plot(:,1)),plot_keep_plot(:,2)];  
    figure(1)
    hold on
    areas_plot(end,2)=1-1E-10;
    semilogx([areas_plot(:,1);areas_plot(end,1)],[1-areas_plot(:,2);1E-10],':','Color',color_mat(dates,:),'Linewidth',3)
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',3,'fontweight','bold','fontsize',16)     
    
    %this is so we can arbitrarily replace our time on approach to scaling
    %and update the initial conditions
    R_zero=min(rads_mean(avg_avg_n_vs_time(:,1)>t_approach_to_scaling));
    R_c_zero=min(rads_mean(avg_avg_n_vs_time(:,1)>t_approach_to_scaling));
    em=b_over_a;
    data_vals(dates,:)=[R_zero,R_c_zero,em];
    time_approch_to_scaling=avg_avg_n_vs_time(avg_avg_n_vs_time(:,1)>t_approach_to_scaling,1);
    %we build our cell array so we can do simulatenous fits
    all_x_data{dates}={rads_mean(avg_avg_n_vs_time(:,1)>t_approach_to_scaling)};
    all_y_data{dates}={time_approch_to_scaling-time_approch_to_scaling(1)};
    all_weights{dates}={1./(time_approch_to_scaling-time_approch_to_scaling(1)+1E-3)};
%     keyboard

    %do not touch anything to do with figure 8
    figure(8)
    plot(avg_a_vs_time_plot(:,1),rads_mean,'o','Linewidth',2,'Color',[103,103,103]/255)
    hold on
    
    times_func_fit_for_rads_data=@(coeffs,rads)...
    (rads-R_zero).*(coeffs(1)*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*rads-coeffs(1)*coeffs(2)*(coeffs(2)+coeffs(1)*coeffs(2)*coeffs(3))*R_zero+2*(coeffs(2)^2)*R_c_zero)/(2*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^2)+...
    coeffs(3)*((coeffs(1)*coeffs(2)*R_zero-coeffs(2)*R_c_zero)^2)*log(((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*rads+coeffs(3)*(coeffs(1)*coeffs(2)*R_zero-coeffs(2)*R_c_zero))/(coeffs(2)*(R_zero-coeffs(3)*R_c_zero)))/((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^3);
    
    coeffs_in=[em,0.05,0.95];
    coeffs_out=nlinfit(rads_mean,avg_avg_n_vs_time(:,1)-avg_avg_n_vs_time(1,1),times_func_fit_for_rads_data,coeffs_in,'Weights',1./(avg_avg_n_vs_time(:,1)-avg_avg_n_vs_time(1,1)+1E-3))      

    
    delta=(1/em)*(em*R_zero-R_c_zero)/(1-em*coeffs_out(2));
    t_0=avg_a_vs_time_plot(1,1);
    
    rads_for_meth_1=rads_mean(avg_a_vs_time_plot(:,1)>t_scaling);
    
    times_for_meth_1=avg_a_vs_time_plot(avg_a_vs_time_plot(:,1)>t_scaling,1);
    R_zero_meth_1=rads_for_meth_1(1);
    t_0_meth_1=times_for_meth_1(1);
    
    plot(times_for_meth_1,rads_for_meth_1,'ok','Linewidth',2)
    method_1_D=@(coeffs,times)(sqrt(coeffs(1)*(times-t_0_meth_1)+R_zero_meth_1^2));
    coeffs_in_1_D=[0.007];
    [coeffs_out_1_D,Rs,Jac,CovB,MSE]=nlinfit(times_for_meth_1,rads_for_meth_1,method_1_D,coeffs_in_1_D);
    coeffs_out_1_D_err=nlparci(coeffs_out_1_D,Rs,Jac);
    [coeffs_out_1_D(1),coeffs_out_1_D-coeffs_out_1_D_err(1,1)]
    
    fit_ys_1_D=method_1_D(coeffs_out_1_D,times_for_meth_1);
    plot(times_for_meth_1,fit_ys_1_D,'--','Linewidth',4,'Color',[0,153/255,0])
    
    method_1=@(coeffs,times)(sqrt(coeffs_out_1_D(1)*(times-t_0)+(R_zero-coeffs(1))^2)+coeffs(1));
    coeffs_in_1=[delta];
    [coeffs_out_1,Rs,Jac,CovB,MSE]=nlinfit(avg_a_vs_time_plot(:,1),rads_mean,method_1,coeffs_in_1);
    %     coeffs_out_1_err=nlparci(coeffs_out_1,Rs,Jac);
    %     [coeffs_out_1(1),coeffs_out_1-coeffs_out_1_err(1,1)]
    
    fit_ys_1=method_1(coeffs_out_1,avg_a_vs_time_plot(:,1));
    plot(avg_a_vs_time_plot(:,1),fit_ys_1,'--','Linewidth',4,'Color',[0,0,204/255])
    
    method_2=@(coeffs,times)(sqrt(coeffs(1)*(times-t_0)+(R_zero-coeffs(2))^2)+coeffs(2));
    coeffs_in_2=[coeffs_out_1_D(1),delta];
    [coeffs_out_2,Rs,Jac,CovB,MSE]=nlinfit(avg_a_vs_time_plot(:,1),rads_mean,method_2,coeffs_in_2);
    %     coeffs_out_2_err=nlparci(coeffs_out_2,Rs,Jac);
    %     [coeffs_out_2(1),coeffs_out_2(1)-coeffs_out_2_err(1,1);...
    %      coeffs_out_2(2),coeffs_out_2(2)-coeffs_out_2_err(2,1)]
    
    fit_ys_2=method_2(coeffs_out_2,avg_a_vs_time_plot(:,1));
    plot(avg_a_vs_time_plot(:,1),fit_ys_2,'--','Linewidth',4,'Color',[204/255,0,0])
    
    
    method_3=@(coeffs,times)(sqrt(coeffs(1)*(times-t_0)+(coeffs(3)-coeffs(2))^2)+coeffs(2));
    coeffs_in_3=[coeffs_out_1_D(1),delta,R_zero];
    [coeffs_out_3,Rs,Jac,CovB,MSE]=nlinfit(avg_a_vs_time_plot(:,1),rads_mean,method_3,coeffs_in_3);
    %     coeffs_out_3_err=nlparci(coeffs_out_3,Rs,Jac);
    %     [coeffs_out_3(1),coeffs_out_3(1)-coeffs_out_3_err(1,1);...
    %      coeffs_out_3(2),coeffs_out_3(2)-coeffs_out_3_err(2,1);...
    %      coeffs_out_3(3),coeffs_out_3(3)-coeffs_out_3_err(3,1)]
    
    fit_ys_3=method_3(coeffs_out_3,avg_a_vs_time_plot(:,1));
    plot(avg_a_vs_time_plot(:,1),fit_ys_3,'--','Linewidth',4,'Color',[204/255,0,102/255])
end

figure(1)
plot_xs=logspace(-4,4);
plot_ys=exp(-0.926*(plot_xs.^1.21));
plot(plot_xs,plot_ys,'--r','Linewidth',4)
plot(plot_xs,exp(-plot_xs),'-.k','Linewidth',3)


%%%Below we do allt he different variety of fits
%the ys and weights are the same for any of the fits, as are out initial
%gueses
yy=vertcat(all_y_data{:});
ww=vertcat(all_weights{:});
coeffs_in=[0.97,0.1,0.97,data_vals(1,1),data_vals(1,2),data_vals(2,1),data_vals(2,2),data_vals(3,1),data_vals(3,2)];
em=coeffs_in(1);

%!!we want to add to the fits for figure 6. Below is the equation as
%written in the paper from Doug this si working as is!! This fits for em, a
%and s simultaneously and hold R_o and R_co constant from the data
times_func_constant_em_fit_a_s=@(coeffs,rads)...
    [(cell2mat(rads{1})-data_vals(1,1)).*(em*coeffs(2)*(coeffs(2)-em*coeffs(2)*coeffs(3))*cell2mat(rads{1})-em*coeffs(2)*(coeffs(2)+em*coeffs(2)*coeffs(3))*data_vals(1,1)+2*(coeffs(2)^2)*data_vals(1,2))/(2*coeffs(2)*(coeffs(2)-em*coeffs(2)*coeffs(3))^2)+...
    coeffs(3)*((em*coeffs(2)*data_vals(1,1)-coeffs(2)*data_vals(1,2))^2)*log(((coeffs(2)-em*coeffs(2)*coeffs(3))*cell2mat(rads{1})+coeffs(3)*(em*coeffs(2)*data_vals(1,1)-coeffs(2)*data_vals(1,2)))/(coeffs(2)*(data_vals(1,1)-coeffs(3)*data_vals(1,2))))/((coeffs(2)-em*coeffs(2)*coeffs(3))^3);...
    (cell2mat(rads{2})-data_vals(2,1)).*(em*coeffs(2)*(coeffs(2)-em*coeffs(2)*coeffs(3))*cell2mat(rads{2})-em*coeffs(2)*(coeffs(2)+em*coeffs(2)*coeffs(3))*data_vals(2,1)+2*(coeffs(2)^2)*data_vals(2,2))/(2*coeffs(2)*(coeffs(2)-em*coeffs(2)*coeffs(3))^2)+...
    coeffs(3)*((em*coeffs(2)*data_vals(2,1)-coeffs(2)*data_vals(2,2))^2)*log(((coeffs(2)-em*coeffs(2)*coeffs(3))*cell2mat(rads{2})+coeffs(3)*(em*coeffs(2)*data_vals(2,1)-coeffs(2)*data_vals(2,2)))/(coeffs(2)*(data_vals(2,1)-coeffs(3)*data_vals(2,2))))/((coeffs(2)-em*coeffs(2)*coeffs(3))^3);...
    (cell2mat(rads{3})-data_vals(3,1)).*(em*coeffs(2)*(coeffs(2)-em*coeffs(2)*coeffs(3))*cell2mat(rads{3})-em*coeffs(2)*(coeffs(2)+em*coeffs(2)*coeffs(3))*data_vals(3,1)+2*(coeffs(2)^2)*data_vals(3,2))/(2*coeffs(2)*(coeffs(2)-em*coeffs(2)*coeffs(3))^2)+...
    coeffs(3)*((em*coeffs(2)*data_vals(3,1)-coeffs(2)*data_vals(3,2))^2)*log(((coeffs(2)-em*coeffs(2)*coeffs(3))*cell2mat(rads{3})+coeffs(3)*(em*coeffs(2)*data_vals(3,1)-coeffs(2)*data_vals(3,2)))/(coeffs(2)*(data_vals(3,1)-coeffs(3)*data_vals(3,2))))/((coeffs(2)-em*coeffs(2)*coeffs(3))^3)];

%this will fit for a and s simultaneously and keep em constant
coeffs_out_constant_em_fit_a_s=nlinfit(all_x_data,cell2mat(yy),times_func_constant_em_fit_a_s,coeffs_in,'Weights',cell2mat(ww));


%!!we want to add to the fits for figure 6. Below is the equation as
%written in the paper from Doug this si working as is!! This fits for em, a
%and s simultaneously and hold R_o and R_co constant from the data
times_func_fit_em_a_s=@(coeffs,rads)...
    [(cell2mat(rads{1})-data_vals(1,1)).*(coeffs(1)*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*cell2mat(rads{1})-coeffs(1)*coeffs(2)*(coeffs(2)+coeffs(1)*coeffs(2)*coeffs(3))*data_vals(1,1)+2*(coeffs(2)^2)*data_vals(1,2))/(2*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^2)+...
    coeffs(3)*((coeffs(1)*coeffs(2)*data_vals(1,1)-coeffs(2)*data_vals(1,2))^2)*log(((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*cell2mat(rads{1})+coeffs(3)*(coeffs(1)*coeffs(2)*data_vals(1,1)-coeffs(2)*data_vals(1,2)))/(coeffs(2)*(data_vals(1,1)-coeffs(3)*data_vals(1,2))))/((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^3);...
    (cell2mat(rads{2})-data_vals(2,1)).*(coeffs(1)*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*cell2mat(rads{2})-coeffs(1)*coeffs(2)*(coeffs(2)+coeffs(1)*coeffs(2)*coeffs(3))*data_vals(2,1)+2*(coeffs(2)^2)*data_vals(2,2))/(2*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^2)+...
    coeffs(3)*((coeffs(1)*coeffs(2)*data_vals(2,1)-coeffs(2)*data_vals(2,2))^2)*log(((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*cell2mat(rads{2})+coeffs(3)*(coeffs(1)*coeffs(2)*data_vals(2,1)-coeffs(2)*data_vals(2,2)))/(coeffs(2)*(data_vals(2,1)-coeffs(3)*data_vals(2,2))))/((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^3);...
    (cell2mat(rads{3})-data_vals(3,1)).*(coeffs(1)*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*cell2mat(rads{3})-coeffs(1)*coeffs(2)*(coeffs(2)+coeffs(1)*coeffs(2)*coeffs(3))*data_vals(3,1)+2*(coeffs(2)^2)*data_vals(3,2))/(2*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^2)+...
    coeffs(3)*((coeffs(1)*coeffs(2)*data_vals(3,1)-coeffs(2)*data_vals(3,2))^2)*log(((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*cell2mat(rads{3})+coeffs(3)*(coeffs(1)*coeffs(2)*data_vals(3,1)-coeffs(2)*data_vals(3,2)))/(coeffs(2)*(data_vals(3,1)-coeffs(3)*data_vals(3,2))))/((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^3)];

%this wll fit for a, s and em simultaneously while using he R_o and R_co
%values found from the data
coeffs_out_fit_em_a_s=nlinfit(all_x_data,cell2mat(yy),times_func_fit_em_a_s,coeffs_in,'Weights',cell2mat(ww));

%%this will fit for s, a and m simultaneously and fit for R_o and R_co for
%%each data set on its own
times_func_fit_all=@(coeffs,rads)...
    [(cell2mat(rads{1})-coeffs(4)).*(coeffs(1)*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*cell2mat(rads{1})-coeffs(1)*coeffs(2)*(coeffs(2)+coeffs(1)*coeffs(2)*coeffs(3))*coeffs(4)+2*(coeffs(2)^2)*coeffs(5))/(2*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^2)+...
    coeffs(3)*((coeffs(1)*coeffs(2)*coeffs(4)-coeffs(2)*coeffs(5))^2)*log(((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*cell2mat(rads{1})+coeffs(3)*(coeffs(1)*coeffs(2)*coeffs(4)-coeffs(2)*coeffs(5)))/(coeffs(2)*(coeffs(4)-coeffs(3)*coeffs(5))))/((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^3);...
    (cell2mat(rads{2})-coeffs(6)).*(coeffs(1)*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*cell2mat(rads{2})-coeffs(1)*coeffs(2)*(coeffs(2)+coeffs(1)*coeffs(2)*coeffs(3))*coeffs(6)+2*(coeffs(2)^2)*coeffs(7))/(2*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^2)+...
    coeffs(3)*((coeffs(1)*coeffs(2)*coeffs(6)-coeffs(2)*coeffs(7))^2)*log(((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*cell2mat(rads{2})+coeffs(3)*(coeffs(1)*coeffs(2)*coeffs(6)-coeffs(2)*coeffs(7)))/(coeffs(2)*(coeffs(6)-coeffs(3)*coeffs(7))))/((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^3);...
    (cell2mat(rads{3})-coeffs(8)).*(coeffs(1)*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*cell2mat(rads{3})-coeffs(1)*coeffs(2)*(coeffs(2)+coeffs(1)*coeffs(2)*coeffs(3))*coeffs(8)+2*(coeffs(2)^2)*coeffs(9))/(2*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^2)+...
    coeffs(3)*((coeffs(1)*coeffs(2)*coeffs(8)-coeffs(2)*coeffs(9))^2)*log(((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*cell2mat(rads{3})+coeffs(3)*(coeffs(1)*coeffs(2)*coeffs(8)-coeffs(2)*coeffs(9)))/(coeffs(2)*(coeffs(8)-coeffs(3)*coeffs(9))))/((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^3)];

%this will fit for s, a and m simultaneously and fit for R_o and R_co for
%each data set on its own
coeffs_out_fit_all=nlinfit(all_x_data,cell2mat(yy),times_func_fit_all,coeffs_in,'Weights',cell2mat(ww));


% data_vals(1,1:2)=coeffs_out_paper(4:5);
% data_vals(2,1:2)=coeffs_out_paper(6:7);
% data_vals(3,1:2)=coeffs_out_paper(8:9);

figure(7)
% %we want to add to the fits for figure 6. Below is the equation as
% %written in the paper from Doug
% times_func_1=@(coeffs,rads)...
%     (rads-data_vals(1,1)).*(coeffs(1)*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*rads-coeffs(1)*coeffs(2)*(coeffs(2)+coeffs(1)*coeffs(2)*coeffs(3))*data_vals(1,1)+2*(coeffs(2)^2)*data_vals(1,2))/(2*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^2)+...
%     coeffs(3)*((coeffs(1)*coeffs(2)*data_vals(1,1)-coeffs(2)*data_vals(1,2))^2)*log(((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*rads+coeffs(3)*(coeffs(1)*coeffs(2)*data_vals(1,1)-coeffs(2)*data_vals(1,2)))/(coeffs(2)*(data_vals(1,1)-coeffs(3)*data_vals(1,2))))/((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^3);
times_func_1=@(coeffs,rads)...
    (rads-coeffs(4)).*(coeffs(1)*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*rads-coeffs(1)*coeffs(2)*(coeffs(2)+coeffs(1)*coeffs(2)*coeffs(3))*coeffs(4)+2*(coeffs(2)^2)*coeffs(5))/(2*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^2)+...
    coeffs(3)*((coeffs(1)*coeffs(2)*coeffs(4)-coeffs(2)*coeffs(5))^2)*log(((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*rads+coeffs(3)*(coeffs(1)*coeffs(2)*coeffs(4)-coeffs(2)*coeffs(5)))/(coeffs(2)*(coeffs(4)-coeffs(3)*coeffs(5))))/((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^3);

x_1=cell2mat(all_x_data{1});
y_1=cell2mat(all_y_data{1});
% w_1=cell2mat(all_weights{1});
% coeffs_out_paper=nlinfit(x_1,y_1,times_func_paper,coeffs_in,'Weights',w_1)
plot(x_1,y_1,symbol_vec{1},'Linewidth',2,'MarkerSize',8,'Color',color_mat(1,:))
hold on
rad_xs_fit=(min(x_1):0.01:10);
fit_ys_constant_em=times_func_1(coeffs_out_constant_em_fit_a_s,rad_xs_fit);
plot(rad_xs_fit,fit_ys_constant_em,'--k','Linewidth',2)
fit_ys_fit_em_a_s=times_func_1(coeffs_out_fit_em_a_s,rad_xs_fit);
plot(rad_xs_fit,fit_ys_fit_em_a_s,'--b','Linewidth',2)
fit_ys_fit_all=times_func_1(coeffs_out_fit_all,rad_xs_fit);
plot(rad_xs_fit,fit_ys_fit_all,'--r','Linewidth',2)

% %we want to add to the fits for figure 6. Below is the equation as
% %written in the paper from Doug
% times_func_2=@(coeffs,rads)...
%     (rads-data_vals(2,1)).*(coeffs(1)*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*rads-coeffs(1)*coeffs(2)*(coeffs(2)+coeffs(1)*coeffs(2)*coeffs(3))*data_vals(2,1)+2*(coeffs(2)^2)*data_vals(2,2))/(2*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^2)+...
%     coeffs(3)*((coeffs(1)*coeffs(2)*data_vals(2,1)-coeffs(2)*data_vals(2,2))^2)*log(((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*rads+coeffs(3)*(coeffs(1)*coeffs(2)*data_vals(2,1)-coeffs(2)*data_vals(2,2)))/(coeffs(2)*(data_vals(2,1)-coeffs(3)*data_vals(2,2))))/((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^3);
times_func_2=@(coeffs,rads)...
    (rads-coeffs(6)).*(coeffs(1)*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*rads-coeffs(1)*coeffs(2)*(coeffs(2)+coeffs(1)*coeffs(2)*coeffs(3))*coeffs(6)+2*(coeffs(2)^2)*coeffs(7))/(2*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^2)+...
    coeffs(3)*((coeffs(1)*coeffs(2)*coeffs(6)-coeffs(2)*coeffs(7))^2)*log(((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*rads+coeffs(3)*(coeffs(1)*coeffs(2)*coeffs(6)-coeffs(2)*coeffs(7)))/(coeffs(2)*(coeffs(6)-coeffs(3)*coeffs(7))))/((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^3);

x_2=cell2mat(all_x_data{2});
y_2=cell2mat(all_y_data{2});
% w_1=cell2mat(all_weights{1});
% coeffs_out_paper=nlinfit(x_1,y_1,times_func_paper,coeffs_in,'Weights',w_1)
plot(x_2,y_2,symbol_vec{2},'Linewidth',2,'MarkerSize',8,'Color',color_mat(2,:))
hold on
rad_xs_fit=(min(x_2):0.01:10);
fit_ys_constant_em=times_func_2(coeffs_out_constant_em_fit_a_s,rad_xs_fit);
plot(rad_xs_fit,fit_ys_constant_em,'--k','Linewidth',2)
fit_ys_fit_em_a_s=times_func_2(coeffs_out_fit_em_a_s,rad_xs_fit);
plot(rad_xs_fit,fit_ys_fit_em_a_s,'--b','Linewidth',2)
fit_ys_fit_all=times_func_2(coeffs_out_fit_all,rad_xs_fit);
plot(rad_xs_fit,fit_ys_fit_all,'--r','Linewidth',2)

% %we want to add to the fits for figure 6. Below is the equation as
% %written in the paper from Doug
% times_func_3=@(coeffs,rads)...
%     (rads-data_vals(3,1)).*(coeffs(1)*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*rads-coeffs(1)*coeffs(2)*(coeffs(2)+coeffs(1)*coeffs(2)*coeffs(3))*data_vals(3,1)+2*(coeffs(2)^2)*data_vals(3,2))/(2*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^2)+...
%     coeffs(3)*((coeffs(1)*coeffs(2)*data_vals(3,1)-coeffs(2)*data_vals(3,2))^2)*log(((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*rads+coeffs(3)*(coeffs(1)*coeffs(2)*data_vals(3,1)-coeffs(2)*data_vals(3,2)))/(coeffs(2)*(data_vals(3,1)-coeffs(3)*data_vals(3,2))))/((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^3);
times_func_3=@(coeffs,rads)...
    (rads-coeffs(8)).*(coeffs(1)*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*rads-coeffs(1)*coeffs(2)*(coeffs(2)+coeffs(1)*coeffs(2)*coeffs(3))*coeffs(8)+2*(coeffs(2)^2)*coeffs(9))/(2*coeffs(2)*(coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^2)+...
    coeffs(3)*((coeffs(1)*coeffs(2)*coeffs(8)-coeffs(2)*coeffs(9))^2)*log(((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))*rads+coeffs(3)*(coeffs(1)*coeffs(2)*coeffs(8)-coeffs(2)*coeffs(9)))/(coeffs(2)*(coeffs(8)-coeffs(3)*coeffs(9))))/((coeffs(2)-coeffs(1)*coeffs(2)*coeffs(3))^3);

x_3=cell2mat(all_x_data{3});
y_3=cell2mat(all_y_data{3});
% w_1=cell2mat(all_weights{1});
% coeffs_out_paper=nlinfit(x_1,y_1,times_func_paper,coeffs_in,'Weights',w_1)
plot(x_3,y_3,symbol_vec{3},'Linewidth',2,'MarkerSize',8,'Color',color_mat(3,:))
hold on
rad_xs_fit=(min(x_3):0.01:10);
fit_ys_constant_em=times_func_3(coeffs_out_constant_em_fit_a_s,rad_xs_fit);
plot(rad_xs_fit,fit_ys_constant_em,'--k','Linewidth',2)
fit_ys_fit_em_a_s=times_func_3(coeffs_out_fit_em_a_s,rad_xs_fit);
plot(rad_xs_fit,fit_ys_fit_em_a_s,'--b','Linewidth',2)
fit_ys_fit_all=times_func_3(coeffs_out_fit_all,rad_xs_fit);
plot(rad_xs_fit,fit_ys_fit_all,'--r','Linewidth',2)

[coeffs_in;coeffs_out_constant_em_fit_a_s;coeffs_out_fit_em_a_s;coeffs_out_fit_all]

keyboard
