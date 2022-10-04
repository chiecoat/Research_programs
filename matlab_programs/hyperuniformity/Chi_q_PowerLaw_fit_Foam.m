% +
% NAME: Chi_q_PowerLaw_fit_Foam
%
% PURPOSE:
%     This program will find the power law scaling for the final decade of
%     the spectral density for foam data.
%
% CATEGORY:
%     Hyperuniformity
%
% CALLING SEQUENCE:
%    Chi_q_PowerLaw_fit_Foam
%
% INPUTS: None at the moment as it is a stand alone program
%
% OPTIONAL INPUTS: (none)
%
% KEYWORD PARAMETERS: (none)
%
% OUTPUTS: Will write a text file with the fit values for the spectral
%          density 

% SIDE EFFECTS: (none)
%
% MODIFICATION HISTORY:
%    written by: A. Chieco, UPenn, November  2018
%-
% clear all
% close all

chiq_path='E:\Chieco\Hyperuniformity\HU foams\fourier space data\foam data';

num={'125','250','500'};
pattern_vec={'points_weighted'};

color_mat=[60,60,60;120,120,120;190,190,190]/255;
pl_color=[0,255,0;0,255,255]/255;
tot_fit=[0,0];
for i1=1:1%numel(pattern_vec)
    pattern=pattern_vec{i1};    
    for i2=1:numel(num)        
        chiq=dlmread([chiq_path '\Chi_q bub_' pattern ' square_imi_' num{i2} ' L_sys_2820.txt']);
        chiq_fit=chiq(2:end,:);
        
        all_keep=chiq_fit(chiq_fit(:,4)>0,:);
        
        x_all=all_keep(:,1);
        x_min=min(x_all);
        x_max=0.6;%10*x_min;
        y_all=all_keep(:,3);
        y_all_weights=real(all_keep(:,4));
        
        
        x_fit=x_all(and(x_all>=x_min,x_all<=x_max));
        y_fit=y_all(and(x_all>=x_min,x_all<=x_max));
        y_weights=1./(y_all_weights(and(x_all>=x_min,x_all<=x_max)).^2);
        
        pl_fit=@(coeffs,x_in)coeffs(1)*x_in.^coeffs(2);
        
        coeffs_in=[10,4];
        [my_fit,Rs,Jac,CovB,MSE]=nlinfit(x_fit,y_fit,pl_fit,coeffs_in,'Weights',y_weights);
        coeffs_err=nlparci(my_fit,Rs,Jac);
        my_err=max(abs([my_fit(2)-coeffs_err(2),my_fit(2)-coeffs_err(2)]));
        pow=my_fit(2);%3.7%
        [my_fit(2),2*my_err/sqrt(numel(y_fit))]
        
        x_plot=logspace(log10(x_min/2),log10(x_max),500);
        y_plot=my_fit(1)*x_plot.^pow;        
        loglog(x_all,y_all,'o-','Color',color_mat(i2,:),'LineWidth',3);
        hold on
        loglog(x_plot,y_plot,'--','Color',pl_color(i1,:),'LineWidth',2);  
        
        y_plot_err_plus=my_fit(1)*x_plot.^(pow+2*my_err/sqrt(numel(y_fit))); 
        loglog(x_plot,y_plot_err_plus,'--','Color',[60 60 60]/255,'LineWidth',2); 
        
        y_plot_err_minus=my_fit(1)*x_plot.^(pow-2*my_err/sqrt(numel(y_fit))); 
        loglog(x_plot,y_plot_err_minus,'--','Color',[60 60 60]/255,'LineWidth',2); 
        
        set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',1)
        xlim([1E-3 3])
        ylim([1E-9 2])
        
        tot_fit=[tot_fit;[x_all,y_all]];
        
    end
    tot_fit=tot_fit(2:end,:);
    [data,index]=sort(tot_fit(:,1));
    tot_fit=tot_fit(index,:);
    
    x_all=tot_fit(:,1);
    x_min=min(x_all);
    x_max=0.6;
    y_all=tot_fit(:,2);
    
    x_fit=x_all(and(x_all>=x_min,x_all<=x_max));
    y_fit=y_all(and(x_all>=x_min,x_all<=x_max));
    y_weights=1./(y_fit.^2);
    
    pl_fit=@(coeffs,x_in)coeffs(1)*x_in.^coeffs(2);
    
    coeffs_in=[10,4.2];
    [my_fit,Rs,Jac,CovB,MSE]=nlinfit(x_fit,y_fit,pl_fit,coeffs_in,'Weights',y_weights);
    coeffs_err=nlparci(my_fit,Rs,Jac);
    my_err=max(abs([my_fit(2)-coeffs_err(2),my_fit(2)-coeffs_err(2)]));
    pow=my_fit(2);
    [my_fit(2),my_err]
    
    x_plot=logspace(log10(x_min/2),log10(x_max),500);
    y_plot=my_fit(1)*x_plot.^pow;
    loglog(x_all,y_all,'o-','Color',color_mat(i2,:),'LineWidth',3);
    hold on
    loglog(x_plot,y_plot,'--','Color',pl_color(i1,:),'LineWidth',2);
    
    y_plot_err_plus=my_fit(1)*x_plot.^(pow+my_err);
    loglog(x_plot,y_plot_err_plus,'--','Color',[60 60 60]/255,'LineWidth',2);
    
    y_plot_err_minus=my_fit(1)*x_plot.^(pow-my_err);
    loglog(x_plot,y_plot_err_minus,'--','Color',[60 60 60]/255,'LineWidth',2);
    
    set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',1)
    xlim([1E-3 3])
    ylim([1E-9 2])
    
    loglog(x_plot,y_plot,'--r','LineWidth',2);
%     filewrite_name=[chiq_path '\power law fits\Chi_power_law bub_' pattern ' square_average L_sys_2820 epsioln=' num2str(pow) ' .txt'];
%     dlmwrite(filewrite_name,[x_plot',y_plot'],'newline','pc')
end
plotedit on