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
loop_start_vec=[1,24,24,40];
loop_end_vec=[69,124,124,200];

color_mat=[0,0,0;...
           204,102,0;...
           0,104,0;...
           204,0,204]/255;

% date_vec={'2021-11-30'};
% 
% num_start=0;
% loop_start_vec=[1];
% loop_end_vec=[60];


for dates=1:numel(date_vec)
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
    p_max=20;
    powers=(1:p_max);
    % n_times=numel(times);
    n_powers=numel(powers);
    
    N_bub_mat=zeros(1,4);
    
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
        
        %we may want to elimiante bubbles based on certain structural
        %characteristics. But for now we will not
        % %       the first point is the interstitial
        %         figure(2)
        %         semilogy(bub_stack(:,5),bub_stack(:,4),'or')
        %         hold on
        %         plot([min_area,min_area],[1,100])
        %         figure(3)
        %         semilogx(bub_stack(:,4),(bub_stack(:,7).^2)./bub_stack(:,4),'or')
        %         hold on
        %         plot([min_area,1E6],[4*pi,4*pi])
        %         plot([min_area,1E6],[5*pi,5*pi])
        %         figure(4)
        %         semilogy(bub_stack(:,6),bub_stack(:,4),'or')
        %         %we go through the process of eliminiating small bubbles with
        %         %different features from the foam. These are most likely 3d or soon
        %         %to be 3d bubbles
        %         elim_rat=find(and(bub_stack(:,4)<1.4,(bub_stack(:,7).^2)./bub_stack(:,4)<5*pi));
        %         elim_ecc=find(and(bub_stack(:,4)<1.7,bub_stack(:,6)>0.9));
        %         elim_other=find(and(and(bub_stack(:,4)<1.4,bub_stack(:,6)>0.7),bub_stack(:,5)<0.87));
        %         elim_vec=sort(unique([elim_rat;elim_ecc;elim_other]));
        %         bub_stack(elim_vec,end)=-1;
        bub_stack=bub_stack(bub_stack(:,end)>0,:);
        
        %plot(bub_stack(:,1),bub_stack(:,2),'ok','Linewidth',2)
        
        n_sides_vec=unique(bub_stack(:,3));
        areas=bub_stack(:,4);
        
        
        area_tot=sum(areas);
        avg_a=mean(areas);
        for i2=1:n_powers
            power=powers(i2);
            avg_a_power=mean(areas.^power);
            A_moment_raw_data(i1-loop_start+1,i2)=avg_a_power/(avg_a^power);
            %now we need to calculate <<n>>_p
            n_data_vec=zeros(1,numel(n_sides_vec));
            for i3=1:numel(n_sides_vec)
                n_sides=n_sides_vec(i3);
                n_sides_data=bub_stack(bub_stack(:,3)==n_sides,:);
                n_data_vec(i3)=(n_sides*sum(n_sides_data(:,4).^power))/(sum(areas.^power));                    
            end
            avg_avg_n_raw_data(i1-loop_start+1,i2)=sum(n_data_vec);
            if i2==1
                F_of_n_mat=[n_sides_vec,transpose(n_data_vec)./n_sides_vec];
            end
        end
        N_bub_mat=[N_bub_mat;[bub_stack(:,1)*0+i1,bub_stack(:,3:4),sort(bub_stack(:,4))/mean(bub_stack(:,4))]];
        
        bub_ids=(0:numel(areas)-1)';
        plot_vals=[sort(bub_stack(:,4)),bub_ids/(numel(areas)-1)];
        plot_keep=plot_vals;
        
        if mod(i1,3)==0
            figure(1)
            semilogx(plot_keep(:,1)/(mean(plot_keep(:,1))),plot_keep(:,2),'Linewidth',1.5,'SeriesIndex',i1)
            set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',1,'fontweight','bold','fontsize',16)
            xlim([0.002 10])
            hold on
            semilogx(plot_keep(:,1)/(mean(plot_keep(:,1))),1-plot_keep(:,2),'Linewidth',1.5,'SeriesIndex',i1)
        end
%         a_cdf=gcf;
%         save_folder='E:\Chieco\Group Meetings\Group Meetings 2021\08-13-2021\area_CDFs';
        %     exportgraphics(a_cdf,[save_folder '\Area_CDF_plot_' num2str(i1) '.png'],'Resolution',200)
        
%         figure(2)
%         bin_edges=(2.5:1:max(bub_stack(:,3))+0.5);
% %         bins=(bin_edges(1:end-1)+bin_edges(2:end))/2;
%         histogram(bub_stack(:,3),bin_edges,'Normalization','pdf','EdgeColor','auto','FaceColor','none','LineWidth',2)
% %         data=histcounts(bub_stack(:,3),bin_edges,'Normalization','pdf');
% %         plot(bins,data,'-o','LineWidth',2,'MarkerSize',6)
%         set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',1,'fontweight','bold','fontsize',16)
%         xlim([2 16])
%         ylim([0 0.5])
%         n_hist=gcf;
%         hold on
%         save_folder='E:\Chieco\Group Meetings\Group Meetings 2021\08-13-2021\n_hists';
%         
%         figure(3)
%         bar(F_of_n_mat(:,1),F_of_n_mat(:,2),1,'LineWidth',2,'EdgeColor','flat','FaceColor','none')
%         hold on
%         set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',1,'fontweight','bold','fontsize',16)
%         xlim([2 16])
%         ylim([0 0.5])
%         hold on
    end
    N_bub_mat=N_bub_mat(2:end,:);
    
    areas_all=N_bub_mat(:,3);
    bub_ids=(0:numel(areas_all)-1)';
    plot_vals=[sort(areas_all),bub_ids/(numel(areas_all)-1)];
    plot_keep=plot_vals;
    
    figure(1)
    semilogx(sort(N_bub_mat(:,end)),plot_keep(:,2),'--k','Linewidth',3)
    hold on
    semilogx(sort(N_bub_mat(:,end)),1-plot_keep(:,2),'--r','Linewidth',3)
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',1,'fontweight','bold','fontsize',16)
    
%     a_cdf=gcf;
%     save_folder='E:\Chieco\Group Meetings\Group Meetings 2021\08-13-2021\area_CDFs';
    % exportgraphics(a_cdf,[save_folder '\Area_CDF_plot_avg.png'],'Resolution',200)
    
    figure(2)
    bin_edges=(0.5:1:max(N_bub_mat(:,2))+0.5);
    histogram(N_bub_mat(:,2),bin_edges,'Normalization','pdf','EdgeColor','k','FaceColor','none','LineWidth',2)
    hold on
    plot([mean(N_bub_mat(:,2)),mean(N_bub_mat(:,2))],[0,0.02],'Color','k')
    xlim([2 16])
    ylim([0 0.5])
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',2,'fontweight','bold','fontsize',16)

    
    figure(3)
    times=unique(N_bub_mat(:,1));
    bin_bars=(2:max(N_bub_mat(:,2)));
    avg_a_vs_time=zeros(numel(times),2);
    avg_a_vs_time(:,1)=5*times;
    avg_avg_n_vs_time=zeros(numel(times),2);
    avg_avg_n_vs_time(:,1)=5*times;
    for t1=1:numel(times)
        F_of_n_vec=bin_bars*0;
        areas_time=N_bub_mat(N_bub_mat(:,1)==t1,:);
        avg_a_vs_time(t1,2)=mean(areas_time(:,3));
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
        if t1==1
            F_to_avg=F_of_n_vec;
        else
            F_to_avg=[F_to_avg;F_of_n_vec];
        end
    end
    bar(bin_bars,mean(F_to_avg),1,'LineWidth',2,'EdgeColor','k','FaceColor','none')
    hold on
    avg_avg_n=sum(bin_bars.*mean(F_to_avg));
    plot([avg_avg_n,avg_avg_n],[0,0.02],'Color','k')
    xlim([2 16])
    ylim([0 0.5])
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',2,'fontweight','bold','fontsize',16)
%     n_hist=gcf;
%     save_folder='E:\Chieco\Group Meetings\Group Meetings 2021\08-13-2021\n_hists';
    % exportgraphics(n_hist,[save_folder '\n_PDF_plot_avg.png'],'Resolution',200)
    % keyboard    
    area_fit_lin=@(coeffs,times) coeffs(1)*(times-coeffs(2));
    coeffs_in=[0.025,100];
    [coeffs_lin,Rs,Jac,CovB,MSE]=nlinfit(avg_a_vs_time(:,1),avg_a_vs_time(:,2),area_fit_lin,coeffs_in);
    
    figure(4)
    plot(avg_a_vs_time(:,1)-coeffs_lin(2),avg_a_vs_time(:,end),'or')
    hold on
    lin_xs=(0:600);
    plot(lin_xs,lin_xs*coeffs_lin(1),'--b','Linewidth',3)
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',2,'fontweight','bold','fontsize',16)
    
    figure(5)
    plot(avg_avg_n_vs_time(:,1)-coeffs_lin(2),avg_avg_n_vs_time(:,end),'ok')
    hold on
    
    mean_identity_vals=zeros(p_max-2,5);
    times=5*((loop_start:num_step:loop_end))'-coeffs_lin(2);
    time_end=2*max(times);
    p_max_plot=5;
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
    
    % figure(6)
    % plot(moment_vec(:,1),moment_vec(:,2),'or','Linewidth',2)
    % hold on
    % plot(moment_vec(:,1),moment_vec(:,3),'ok','Linewidth',2)
    
    figure(9)
    x_vals=(1:p_max);
    x_plot_vals=(1.1:0.1:30);
    y_plot_vals=(x_plot_vals.^2)./(x_plot_vals.^2-1);
    loglog(x_plot_vals,y_plot_vals-1,'linewidth',2,'Color',[153 153 153]/255)
    hold on
    A_moment_vals=(moment_vec(2:end-1,2).^2)./(moment_vec(1:end-2,2).*moment_vec(3:end,2));
    n_moment_vals=(moment_vec(2:end-1,3)-6)./(moment_vec(1:end-2,3)-6);
    y_vals=n_moment_vals.*A_moment_vals;
    error=std(y_vals)/sqrt(numel(y_vals));
    loglog(mean_identity_vals(:,1),mean_identity_vals(:,2)-1,'o','Linewidth',2,'MarkerSize',8,'Color',color_mat(dates,:))
    % plot(x_vals(2:end-1),y_vals,'or','Linewidth',2);
    % plot(mean_identity_vals(:,1),mean_identity_vals(:,3),'og','Linewidth',2)
    % plot(mean_identity_vals(:,1),mean_identity_vals(:,4),'ob','Linewidth',2)    
%     ylim([0.95 1.35])
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.03],'LineWidth',3,'FontSize',10)
        set(gcf,'position',[0,0,1430,1000])
    
%     close Figure 1
    close Figure 2
    close Figure 3
    close Figure 4
%     keyboard
end

