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

% date_vec={'2021-11-30'};
% num_start=0;
% loop_start_vec=[1];
% loop_end_vec=[60];
symbol_vec={'o','+','x'};
file_name_ext={'C','B','A'};
m_size=[8,12,12];

k_0_fit=[0.0265,0.0237,0.0237];


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
        for i2=1:n_powers
            power=powers(i2);
            avg_a_power=mean(areas.^power);
            A_moment_raw_data(i1-loop_start+1,i2)=avg_a_power/(avg_a^power);
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
        N_bub_mat=[N_bub_mat;[bub_stack(:,1)*0+i1,bub_stack(:,3:4),bub_stack(:,4)/2.1]];
        
        
        %confirmed 13 sided bubble in dates=1, i1=31
        %confirmed 14 sided bubble in dates=1, i1=34
        
        bub_ids=(0:numel(areas)-1)';
        plot_vals=[sort(bub_stack(:,4)),bub_ids/(numel(areas)-1)];
        plot_keep=plot_vals;
        
        if mod(i1,3)==0
            %adding individual curves for the area distribution
            figure(1) 
            semilogy(plot_keep(:,1)/(mean(plot_keep(:,1))),1-plot_keep(:,2),'LineWidth',1,'Color',[190 190 190]/255)
            set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',3,'fontweight','bold','fontsize',16)
            xlim([0 8])
            ylim([2E-4 1])
            hold on            
%             %adding individual curves for the side number distribution
%             figure(2)
%             loop_bin_edges=(0.5:1:max(bub_stack(:,3))+0.5);
%             histogram(bub_stack(:,3),loop_bin_edges,'Normalization','pdf','EdgeColor',[153,153,153]/255,'FaceColor','none','LineWidth',1,'DisplayStyle','stairs')
%             hold on
%             %adding individual curves for F(n) distribution
%             figure(3)
%             x_vals=F_of_n_mat(:,1);
%             y_vals=F_of_n_mat(:,2);
%             vals_plot=zeros(2*numel(y_vals),2);
%             for y1=1:2:numel(vals_plot(:,1))-1
%                 counter=(y1-1)/2+1;
%                 vals_plot(y1:y1+1,2)=y_vals(counter);
%                 vals_plot(y1,1)=x_vals(counter)-0.5;
%                 vals_plot(y1+1,1)=x_vals(counter)+0.5;
%             end
%             plot(vals_plot(:,1),vals_plot(:,2),'Color',[153,153,153]/255,'Linewidth',1)
%             hold on
             N_bub_mat_plot=[N_bub_mat_plot;[bub_stack(:,1)*0+i1,bub_stack(:,3:4),sort(bub_stack(:,4))/mean(bub_stack(:,4))]];
        end
%         [dates,numel(areas)]
%         [2*sqrt(mean(areas)/pi)*mm_conv,sqrt(mean(areas))*mm_conv]
    end
%     continue
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
    
    
    %plot histogram for F(n)
    figure(3)
    times=5*unique(N_bub_mat(:,1));
    bin_bars=(2:20);
    %this is to find the average area
    avg_a_vs_time=zeros(numel(times),2);
    avg_a_vs_time(:,1)=times;
    %this is to find the average firs moment
    avg_asq_over_avg_a_sq_vs_time=zeros(numel(times),2);
    avg_asq_over_avg_a_sq_vs_time(:,1)=times;
    %this is to find the average sides number
    avg_n_vs_time=zeros(numel(times),3);
    avg_n_vs_time(:,1)=times;
    %this is to find the average area weighted side number
    avg_avg_n_vs_time=zeros(numel(times),2);
    avg_avg_n_vs_time(:,1)=times;
    data_to_write=[0,0,0,0];
    min_data_to_plot=zeros(numel(times),4);
    for t1=1:numel(times)
        F_of_n_vec=bin_bars*0;
        areas_time=N_bub_mat(5*N_bub_mat(:,1)==times(t1),:);
        avg_a_vs_time(t1,2:3)=[mean(areas_time(:,3)),mean(areas_time(:,3))/sqrt(numel(areas_time(:,3)))];
        if or(t1==1, t1==numel(times))==1
            numel(areas_time(:,3))            
        end
        avg_asq_over_avg_a_sq_vs_time(t1,2)=mean(areas_time(:,3).^2)/(mean(areas_time(:,3))^2);
        avg_n_vs_time(t1,2)=mean(areas_time(:,2));        
        area_tot=sum(areas_time(:,3));
        min_data_check=areas_time(areas_time(:,4)<1,:);
        min_data_to_plot(t1,:)=[times(t1),min(areas_time(:,3)),numel(min_data_check(:,1))/numel(areas_time(:,4)),mean(min_data_check(:,2))];
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
        if t1==1
            F_to_avg=F_of_n_vec;
        else
            F_to_avg=[F_to_avg;F_of_n_vec];
        end
        data_to_write=[data_to_write;[times(t1)+areas_time(:,3)*0,areas_time(:,3),areas_time(:,3)/mean(areas_time(:,3)),areas_time(:,2)]];
    end
    y_vals=mean(F_to_avg);
    y_err=std(F_to_avg)/sqrt(numel(times));
    x_vals=bin_bars;
    vals_plot=zeros(2*numel(y_vals),2);
    for y1=1:2:numel(vals_plot(:,1))-1
        counter=(y1-1)/2+1;
        vals_plot(y1:y1+1,2)=y_vals(counter);
        vals_plot(y1,1)=bin_bars(counter)-0.5;
        vals_plot(y1+1,1)=bin_bars(counter)+0.5;
    end
    
    
    data=zeros(19,6);
    n_bub_tot=numel(N_bub_mat(:,1))
    a_tot_check=sum(N_bub_mat(:,3));
    n_photos=max(N_bub_mat(:,1))-min(N_bub_mat(:,1));
    for i1_check=2:20
        n_bubs=N_bub_mat(N_bub_mat(:,2)==i1_check,:);
        num_bubs=numel(n_bubs(:,1));
        data(i1_check-2+1,1:2)=[i1_check,numel(n_bubs(:,1))];
        if num_bubs==0
            continue
        else
            p_of_n=num_bubs/n_bub_tot;
            data(i1_check-2+1,3:4)=[p_of_n,p_of_n/sqrt(num_bubs)];
        end
    end
    data(:,5:6)=[y_vals',y_err'];
    
    
    if dates==1
        n_side_mat=[N_bub_mat(:,2),zeros(numel(N_bub_mat(:,2)),1)+dates];
        F_of_n_vals_total=[vals_plot,zeros(numel(vals_plot(:,1)),1)+dates];
        avg_avg_n_vec=sum(bin_bars.*mean(F_to_avg));
        data_tot=[data,zeros(numel(data(:,1)),1)+dates];
    else
        n_side_mat=[n_side_mat;N_bub_mat(:,2),zeros(numel(N_bub_mat(:,2)),1)+dates];
        F_of_n_vals_total=[F_of_n_vals_total;[vals_plot,zeros(numel(vals_plot(:,1)),1)+dates]];
        avg_avg_n_vec=[avg_avg_n_vec,sum(bin_bars.*mean(F_to_avg));]
        data_tot=[data_tot;data,zeros(numel(data(:,1)),1)+dates];
    end
    
    avg_a_vs_time_plot=avg_a_vs_time(avg_a_vs_time(:,2)>0,:);
    area_fit_lin=@(coeffs,times) coeffs(1)*(times-coeffs(2));
    coeffs_in=[0.025,100];
    [coeffs_lin,Rs,Jac,CovB,MSE]=nlinfit(avg_a_vs_time_plot(:,1),avg_a_vs_time_plot(:,2),area_fit_lin,coeffs_in);
    
    figure(4)
    if dates==1
        errorbar(avg_a_vs_time_plot(:,1)-coeffs_lin(2),avg_a_vs_time_plot(:,2),avg_a_vs_time_plot(:,3),symbol_vec{dates},'Linewidth',2,'MarkerSize',8,'Color',color_mat(dates,:),'capSize',0)
        [dates,avg_a_vs_time_plot(1,2),avg_a_vs_time_plot(end,2),numel(avg_a_vs_time_plot(:,2))]
        keyboard
    else
        plot(avg_a_vs_time_plot(:,1)-coeffs_lin(2),avg_a_vs_time_plot(:,2),symbol_vec{dates},'Linewidth',2,'MarkerSize',8,'Color',color_mat(dates,:))
        [dates,avg_a_vs_time_plot(1,2),avg_a_vs_time_plot(end,2),numel(avg_a_vs_time_plot(:,2))]
        keyboard
    end
     hold on
    lin_xs=(0:5:1E3);
%     slope=2*k_0_fit(dates)*(1/mean(avg_asq_over_avg_a_sq_vs_time(:,end)))*(mean(avg_avg_n_vs_time(:,end))-6)
%     plot(lin_xs,lin_xs*slope,'-.k','Linewidth',3)
    plot(lin_xs,lin_xs*coeffs_lin(1),'-.k','Linewidth',3)
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',2,'fontweight','bold','fontsize',16)
    coeffs_lin_err=nlparci(coeffs_lin,Rs,Jac);
    [coeffs_lin(1),coeffs_lin(1)-coeffs_lin_err(1,1)]
    
    figure(11)
    plot(min_data_to_plot(:,1)-coeffs_lin(2),min_data_to_plot(:,3),symbol_vec{dates},'Linewidth',2,'MarkerSize',8,'Color',color_mat(dates,:))
    hold on
    plot(min_data_to_plot(:,1)-coeffs_lin(2),min_data_to_plot(:,4),symbol_vec{dates},'Linewidth',2,'MarkerSize',8,'Color',color_mat(dates,:))
    
    
    data_to_write=data_to_write(2:end,:);
    data_to_write(:,1)=data_to_write(:,1)-coeffs_lin(2);    
    table_to_write = array2table(data_to_write);
    table_to_write.Properties.VariableNames(1:4) = {'time(mins)','Area(mm^2)','A/<A>','n'};
    writetable(table_to_write,['E:\Chieco\My Papers\foam moments\Foam_' file_name_ext{dates} '_data.csv'])
    
    figure(5)
    plot(avg_avg_n_vs_time(:,1)-coeffs_lin(2),avg_avg_n_vs_time(:,end),symbol_vec{dates},'Linewidth',2,'MarkerSize',8,'Color',color_mat(dates,:))
    hold on
    plot([0,1E3],[mean(avg_avg_n_vs_time(:,end)),mean(avg_avg_n_vs_time(:,end))],':','Linewidth',2,'Color',color_mat(dates,:));
    plot(avg_asq_over_avg_a_sq_vs_time(:,1)-coeffs_lin(2),avg_asq_over_avg_a_sq_vs_time(:,2),symbol_vec{dates},'Linewidth',2,'MarkerSize',8,'Color',color_mat(dates,:))
    plot([0,1E3],[mean(avg_asq_over_avg_a_sq_vs_time(:,end)),mean(avg_asq_over_avg_a_sq_vs_time(:,end))],'-.','Linewidth',2,'Color',color_mat(dates,:));
    plot(avg_n_vs_time(:,1)-coeffs_lin(2),avg_n_vs_time(:,end),symbol_vec{dates},'Linewidth',2,'MarkerSize',8,'Color',color_mat(dates,:))
    plot([0,1E3],[mean(avg_n_vs_time(:,end)),mean(avg_n_vs_time(:,end))],'--','Linewidth',2,'Color',color_mat(dates,:));
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',2,'fontweight','bold','fontsize',16)
    
    k_0=coeffs_lin(1)*mean(avg_asq_over_avg_a_sq_vs_time(:,end))/(2*(mean(avg_avg_n_vs_time(:,end))-6))
    area_moment=mean(avg_asq_over_avg_a_sq_vs_time(:,end))
    area_moment_err=std(avg_asq_over_avg_a_sq_vs_time(:,end))
    avg_avg_n=mean(avg_avg_n_vs_time(:,end))
    avg_avg_n_err=std(avg_avg_n_vs_time(:,end))
%     keyboard
    
    if dates==1
        mean_mat=[mean(avg_avg_n_vs_time(:,end)),std(avg_avg_n_vs_time(:,end))/sqrt(numel(avg_avg_n_vs_time(:,end))),...
                  mean(avg_n_vs_time(:,end)),std(avg_n_vs_time(:,end))/sqrt(numel(avg_n_vs_time(:,end))),...
                  mean(avg_asq_over_avg_a_sq_vs_time(:,end)),std(avg_asq_over_avg_a_sq_vs_time(:,end))/sqrt(numel(mean(avg_asq_over_avg_a_sq_vs_time(:,end))))];
    else
        mean_mat=[mean_mat;...
                  mean(avg_avg_n_vs_time(:,end)),std(avg_avg_n_vs_time(:,end))/sqrt(numel(avg_avg_n_vs_time(:,end))),...
                  mean(avg_n_vs_time(:,end)),std(avg_n_vs_time(:,end))/sqrt(numel(avg_n_vs_time(:,end))),...
                  mean(avg_asq_over_avg_a_sq_vs_time(:,end)),std(avg_asq_over_avg_a_sq_vs_time(:,end))/sqrt(numel(mean(avg_asq_over_avg_a_sq_vs_time(:,end))))]; 
    end
    
    mean_identity_vals=zeros(p_max-2,5);
    times=5*((loop_start:num_step:loop_end))'-coeffs_lin(2);
    time_end=2*max(times);
    p_max_plot=5;
    vals_to_plot_vs_p=[0,0,0];
    for i3=1:p_max
        figure(6)
        subplot(1,2,1)
        semilogy(times,A_moment_raw_data(:,i3),'-o','SeriesIndex',i3,'Linewidth',2)
        hold on
        plot([0,time_end],[mean(A_moment_raw_data(:,i3)),mean(A_moment_raw_data(:,i3))],'--k','Linewidth',2)
        xlim([0 400])
        subplot(1,2,2)
        plot(times,avg_avg_n_raw_data(:,i3),'-o','SeriesIndex',i3,'Linewidth',2)
        hold on
        plot([0,time_end],[mean(avg_avg_n_raw_data(:,i3)),mean(avg_avg_n_raw_data(:,i3))],'--k','Linewidth',2)
        xlim([0 400])
        if and(i3>1,i3<p_max)==1
            A_moment_vals=(A_moment_raw_data(:,i3).^2)./(A_moment_raw_data(:,i3-1).*A_moment_raw_data(:,i3+1));
            n_moment_vals=(avg_avg_n_raw_data(:,i3)-6)./(avg_avg_n_raw_data(:,i3-1)-6);
            y_vals=n_moment_vals.*A_moment_vals;
            y_errs=mean(A_moment_vals)*mean(n_moment_vals)*sqrt((std(A_moment_vals)/mean(A_moment_vals))^2+(std(n_moment_vals)/mean(n_moment_vals))^2);
            mean_identity_vals(i3-1,:)=[i3,mean(y_vals),mean(y_vals(and(times>50,times<250))),mean(y_vals(and(times>100,times<276))),y_errs]; 
            vals_to_plot_vs_p=[i3,mean(A_moment_raw_data(:,i3)./(avg_a_vs_time_plot.^i3)),mean(avg_avg_n_raw_data(:,i3))];
            if and(i3>1,i3<=p_max_plot)==1
                figure(7)
                subplot(1,2,1)
                semilogy(times,A_moment_vals,'-o','SeriesIndex',i3,'Linewidth',2)
                set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.03],'LineWidth',3,'FontSize',10)
%                 semilogy(times,A_moment_vals,'-o','Linewidth',2,'Color',color_mat(dates,:))
                hold on
                plot([0,time_end],[mean(A_moment_vals),mean(A_moment_vals)],'--','SeriesIndex',i3,'Linewidth',2)
                xlim([0 400])
                subplot(1,2,2)
                plot(times,n_moment_vals,'-o','SeriesIndex',i3,'Linewidth',2)
                set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.03],'LineWidth',3,'FontSize',10)
%                 plot(times,n_moment_vals,'-o','Color',color_mat(dates,:),'Linewidth',2)
                hold on
                plot([0,time_end],[mean(n_moment_vals),mean(n_moment_vals)],'--','SeriesIndex',i3,'Linewidth',2)
                xlim([0 400])
                figure(8)
                loglog(times,y_vals,'-o','Linewidth',1,'Color',color_mat(dates,:));
                hold on
%                 plot(times,times*0+mean(y_vals),'--','SeriesIndex',i3,'Linewidth',2)
            end
        end
    end
    
    
    n_sides_vec=unique(N_bub_mat(:,2));
    area_tot=sum(areas_all);
    avg_a=mean(areas_all);
    
    moment_vec=zeros(p_max,3);
    for i3=1:p_max
        power=i3;
        moment_vec(i3,1)=power;
        avg_a_power=mean(areas_all.^power);
        moment_vec(i3,2)=avg_a_power;%/(avg_a^power);
        %now we need to calculate <<n>>_p
        n_data_vec=zeros(1,numel(n_sides_vec));
        for i4=1:numel(n_sides_vec)
            n_sides=n_sides_vec(i4);
            n_sides_data=N_bub_mat(N_bub_mat(:,2)==n_sides,:);
            n_data_vec(i4)=(n_sides*sum(n_sides_data(:,3).^power))/(sum(areas_all.^power));
        end
        moment_vec(i3,3)=sum(n_data_vec);
    end
    
    %we plot the identity
    figure(9)
    A_moment_vals=(moment_vec(2:end-1,2).^2)./(moment_vec(1:end-2,2).*moment_vec(3:end,2));
    n_moment_vals=(moment_vec(2:end-1,3)-6)./(moment_vec(1:end-2,3)-6);
    y_vals=n_moment_vals.*A_moment_vals;
    error=std(y_vals)/sqrt(numel(y_vals));
    loglog(mean_identity_vals(:,1),mean_identity_vals(:,2)-1,symbol_vec{dates},'Linewidth',2,'MarkerSize',m_size(dates),'Color',color_mat(dates,:))
    hold on
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.03],'LineWidth',3,'FontSize',10)
    set(gcf,'position',[0,0,1430,1000])
    
    figure(10)
    %this is <A^p>/<A>^p vs p
    loglog(vals_to_plot_vs_p(:,1),vals_to_plot_vs_p(:,2),symbol_vec{dates},'Linewidth',2,'MarkerSize',m_size(dates),'Color',color_mat(dates,:))
    hold on
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.03],'LineWidth',3,'FontSize',10)
    
%     figure(11)
%     %this is <<n>>_p vs p
%     loglog(vals_to_plot_vs_p(:,1),vals_to_plot_vs_p(:,3),symbol_vec{dates},'Linewidth',2,'MarkerSize',m_size(dates),'Color',color_mat(dates,:))
%     hold on
%     set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.03],'LineWidth',3,'FontSize',10)
        
    
    rads_keep=sqrt(plot_keep_plot(:,1));
	if dates==1
        areas_plot_keep_total=[plot_keep_plot(:,1)/mean(plot_keep_plot(:,1)),plot_keep_plot(:,2),plot_keep_plot(:,2)*0+dates];
        rads_plot_keep_total=[rads_keep/mean(rads_keep),plot_keep_plot(:,2),plot_keep_plot(:,2)*0+dates];
    else
        areas_plot_keep_total=[areas_plot_keep_total;...
                               plot_keep_plot(:,1)/mean(plot_keep_plot(:,1)),plot_keep_plot(:,2),plot_keep_plot(:,2)*0+dates];
        rads_plot_keep_total=[rads_plot_keep_total;...
                              rads_keep/mean(rads_keep),plot_keep_plot(:,2),plot_keep_plot(:,2)*0+dates];
    end
%     close Figure 2
    close Figure 3
    close Figure 6
    close Figure 7
    close Figure 8
    %keyboard
end

for dates_plot=1:3
    areas_plot=areas_plot_keep_total(areas_plot_keep_total(:,3)==dates_plot,:);
    rads_plot=rads_plot_keep_total(rads_plot_keep_total(:,3)==dates_plot,:);
    data_avg_keep=data_tot(data_tot(:,end)==dates_plot,:);
     figure(1)
%     semilogx(areas_plot(:,1),areas_plot(:,2),'--','Color',color_mat(dates_plot,:),'Linewidth',3)
    hold on
    F1_mat=[areas_plot(:,1),1-areas_plot(:,2)];
    if dates_plot==1
        dlmwrite('E:\Chieco\Presentations\March Meeting 2022\area cdf\area_cdf_C.txt',F1_mat,'delimiter',',')
    elseif dates_plot==2
        dlmwrite('E:\Chieco\Presentations\March Meeting 2022\area cdf\area_cdf_B.txt',F1_mat,'delimiter',',')
    elseif dates_plot==3
        dlmwrite('E:\Chieco\Presentations\March Meeting 2022\area cdf\area_cdf_A.txt',F1_mat,'delimiter',',')
    end
    areas_plot(end,2)=1-1E-10;
    semilogx(areas_plot(:,1),1-areas_plot(:,2),':','Color',color_mat(dates_plot,:),'Linewidth',3)
    semilogx(areas_plot(end-5:end-1,1),1-areas_plot(end-5:end-1,2),symbol_vec{dates_plot},'Color',color_mat(dates_plot,:),'Linewidth',2)
%     semilogx(2.1/[],1-areas_plot(:,2),':','Color',color_mat(dates_plot,:),'Linewidth',3)
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',3,'fontweight','bold','fontsize',16)
    %this is for the side number distributions
    figure(2)
    n_side_keep=n_side_mat(n_side_mat(:,end)==dates_plot,1);
    bin_edges=(0.5:1:20.5);
    data_to_avg=histogram(n_side_keep,bin_edges,'Normalization','pdf','EdgeColor',color_mat(dates_plot,:),'FaceColor','none','LineWidth',3,'DisplayStyle','stairs');
    hold on
    if dates_plot==1
        n_data_avg_mat=zeros(numel(bin_edges)-1,4);
        n_data_avg_mat(:,1)=(bin_edges(:,1:end-1)+bin_edges(:,2:end))/2;
    end
    n_data_avg_mat(:,dates_plot+1)=transpose(data_to_avg.Values);
    errorbar(data_avg_keep(:,1),data_avg_keep(:,3),data_avg_keep(:,4),'.','Linewidth',3,'Capsize',0,'Color',color_mat(dates_plot,:))
    %this is for the F(n) plots
    figure(3)
    F_of_n_vals_keep=F_of_n_vals_total(F_of_n_vals_total(:,3)==dates_plot,1:2);
    plot(F_of_n_vals_keep(:,1),F_of_n_vals_keep(:,2),'Linewidth',3,'Color',color_mat(dates_plot,:))
    hold on
    errorbar(data_avg_keep(:,1),data_avg_keep(:,5),data_avg_keep(:,6),'.','Linewidth',3,'Capsize',0,'Color',color_mat(dates_plot,:))
    figure(6)
    %this is the PDF of the normalized radii
    rad_data=histogram(rads_plot(:,1),'Normalization','PDF','EdgeColor',color_mat(dates_plot,:),'FaceColor','none','LineWidth',3,'DisplayStyle','stairs');
    hold on
    rad_bins=rad_data.BinEdges;
    x_data=(rad_bins(1:end-1)+rad_bins(2:end))/2;
    y_data=rad_data.Values;
    if dates_plot==1
        data_lognorm_tot=[x_data',y_data'];
    else
        data_lognorm_tot=[data_lognorm_tot;x_data',y_data'];
    end
%     log_normal_coeffs=mle(rads_plot(:,1),'Distribution','Lognormal');
%     plot(x_data,lognpdf(x_data,log_normal_coeffs(1),log_normal_coeffs(2)),'g')
end

figure(1)
plot_xs=logspace(-4,4);
plot_ys=exp(-0.926*(plot_xs.^1.21));
plot(plot_xs,plot_ys,'--r','Linewidth',4)
plot(plot_xs,exp(-plot_xs),'-.k','Linewidth',3)

figure(2)
%we average the total values
bin_edges=(0.5:1:20.5);
data_to_avg=histogram(n_side_mat(:,1),bin_edges,'Normalization','pdf','EdgeColor','k','FaceColor','none','LineWidth',3,'DisplayStyle','stairs');
  
for N1=1:numel(n_data_avg_mat(:,1))
    if N1==1
        n_avg_plot_mat=[n_data_avg_mat(N1,1),mean(n_data_avg_mat(N1,2:4)),(max(n_data_avg_mat(N1,2:4))-min(n_data_avg_mat(N1,2:4)))/2];
        n_bar_plot_mat=[n_data_avg_mat(N1,1)-0.5,mean(n_data_avg_mat(N1,2:4));...
                        n_data_avg_mat(N1,1)+0.5,mean(n_data_avg_mat(N1,2:4))];
    else
        n_avg_plot_mat=[n_avg_plot_mat;
                        n_data_avg_mat(N1,1),mean(n_data_avg_mat(N1,2:4)),(max(n_data_avg_mat(N1,2:4))-min(n_data_avg_mat(N1,2:4)))/2];
        n_bar_plot_mat=[n_bar_plot_mat;...
                        n_data_avg_mat(N1,1)-0.5,mean(n_data_avg_mat(N1,2:4));...
                        n_data_avg_mat(N1,1)+0.5,mean(n_data_avg_mat(N1,2:4))];
    end
    
end
plot(n_bar_plot_mat(:,1),n_bar_plot_mat(:,2),'k','Linewidth',3)
errorbar(n_avg_plot_mat(:,1),n_avg_plot_mat(:,2),n_avg_plot_mat(:,3),'.k','Linewidth',3,'Capsize',0)
plot([mean(n_side_mat(:,1)),mean(n_side_mat(:,1))],[0,0.02],'Color','k')
xlim([2 16])
ylim([0 0.5])
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',3,'fontweight','bold','fontsize',16)

adam_data_n=[3,0.013,0.001;
           4,0.083,0.002;
           5,0.304,0.005;
           6,0.316,0.005;
           7,0.178,0.003;
           8,0.076,0.002;
           9,0.022,0.001;
           10,0.006,0.0006;
           11,0.0005,0.0002];

errorbar(adam_data_n(:,1),adam_data_n(:,2),adam_data_n(:,3),'or','CapSize',0) 

figure(3)

xlim([2 16])
ylim([0 0.5])
F_of_n_avg=[F_of_n_vals_total(1:38,1:2),F_of_n_vals_total(39:76,2),F_of_n_vals_total(77:114,2)];
F_avg_plot=F_of_n_avg(:,1:3);
F_avg_plot(:,2:3)=F_avg_plot(:,2:3)*0;
for F1=1:numel(F_of_n_avg(:,1))
    F_avg_plot(F1,2)=mean(F_of_n_avg(F1,2:4));
end
plot(F_avg_plot(:,1),F_avg_plot(:,2),'k','Linewidth',3)
hold on
%now we need to make the figure have error bats
for F2=1:2:numel(F_of_n_avg(:,1))-1
    if F2==1
        F_of_n_err=[F_of_n_avg(F2)+0.5,mean(F_of_n_avg(F2,2:4)),(max(F_of_n_avg(F2,2:4))-min(F_of_n_avg(F2,2:4)))/2];
    else
        F_of_n_err=[F_of_n_err;F_of_n_avg(F2)+0.5,mean(F_of_n_avg(F2,2:4)),(max(F_of_n_avg(F2,2:4))-min(F_of_n_avg(F2,2:4)))/2];
    end      
end
errorbar(F_of_n_err(:,1),F_of_n_err(:,2),F_of_n_err(:,3),'.k','Linewidth',3,'Capsize',0)
plot([mean(avg_avg_n_vec),mean(avg_avg_n_vec)],[0,0.02],'Color','k')

adam_data=[3,0.0009,0.00006;
           4,0.034,0.006;
           5,0.173,0.007;
           6,0.326,0.009;
           7,0.259,0.007;
           8,0.141,0.006;
           9,0.049,0.004;
           10,0.016,0.002;
           11,0.001,0.0005];
errorbar(adam_data(:,1),adam_data(:,2),adam_data(:,3),'or','CapSize',0)           
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',3,'fontweight','bold','fontsize',16) 


figure(5)
F_of_n_mean_weights=1./(mean_mat(:,2).^2);
F_of_n_mean_error=sqrt(1/sum(F_of_n_mean_weights));
F_of_n_mean=sum(mean_mat(:,1).*F_of_n_mean_weights)/sum(F_of_n_mean_weights);
plot([0,1E3],[F_of_n_mean,F_of_n_mean],':','Linewidth',2,'Color','k')
[F_of_n_mean,sqrt(3)*F_of_n_mean_error]

n_mean_weights=1./(mean_mat(:,4).^2);
n_mean_error=sqrt(1/sum(n_mean_weights));
n_mean=sum(mean_mat(:,3).*n_mean_weights)/sum(n_mean_weights);
plot([0,1E3],[n_mean,n_mean],'--','Linewidth',2,'Color','k');
[n_mean,sqrt(3)*n_mean_error]

Asq_o_A_mean_weights=1./(mean_mat(:,6).^2);
Asq_o_A_mean_error=sqrt(1/sum(Asq_o_A_mean_weights));
Asq_o_A_mean=sum(mean_mat(:,5).*Asq_o_A_mean_weights)/sum(Asq_o_A_mean_weights);
plot([0,1E3],[Asq_o_A_mean,Asq_o_A_mean],'-.','Linewidth',2,'Color','k');
[Asq_o_A_mean,sqrt(3)*Asq_o_A_mean_error]

figure(6)
%now we fit the data to a log normal distribution
[logn_sort,ids]=sort(data_lognorm_tot(:,1));
log_norm_fit=data_lognorm_tot(ids,:);


log_normal_fit=@(sigma,xs) (1./(sigma(1)*xs*sqrt(2*pi))).*exp(-((log(xs)+sigma(1)^2/2).^2)/(2*sigma(1)^2));
coeffs_in_lognormal=[0.5];
[sigma_out,Rs,Jac,CovB,MSE]=nlinfit(log_norm_fit(:,1)',log_norm_fit(:,2)',log_normal_fit,coeffs_in_lognormal,'Weights',log_norm_fit(:,2)');
x_plot_data=(0:max(log_norm_fit(:,1))/100:2*max(log_norm_fit(:,1)));
y_data_fit=(1./(sigma_out*x_plot_data*sqrt(2*pi))).*exp(-((log(x_plot_data)+sigma_out^2/2).^2)/(2*sigma_out^2));
plot(x_plot_data,y_data_fit,'--k','Linewidth',4)

weibull_coeffs=mle(sort(rads_plot_keep_total(:,1)),'Distribution','wbl');
plot(x_plot_data,wblpdf(x_plot_data,weibull_coeffs(1),weibull_coeffs(2)),'r','Linewidth',2)
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',3,'fontweight','bold','fontsize',16)

log_normal_coeffs=mle(sort(rads_plot_keep_total(:,1)),'Distribution','Lognormal');
plot(x_plot_data,lognpdf(x_plot_data,log_normal_coeffs(1),log_normal_coeffs(2)),'g')
sigma_out=log_normal_coeffs(2);
y_data_fit=(1./(sigma_out*x_plot_data*sqrt(2*pi))).*exp(-((log(x_plot_data)+sigma_out^2/2).^2)/(2*sigma_out^2));
plot(x_plot_data,y_data_fit,'--k','Linewidth',4)


figure(9)
x_plot_vals=(1.1:0.1:30);
y_plot_vals=(x_plot_vals.^2)./(x_plot_vals.^2-1);
loglog(x_plot_vals,y_plot_vals-1,'linewidth',2,'Color','k')
