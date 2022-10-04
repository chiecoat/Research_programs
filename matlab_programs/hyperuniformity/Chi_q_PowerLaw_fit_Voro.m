% +
% NAME: Chi_q_PowerLaw_fit_Voro
%
% PURPOSE:
%     This program will find the power law scaling for the final decade of
%     the spectral density for the voronoi patterns.
%
% CATEGORY:
%     Hyperuniformity
%
% CALLING SEQUENCE:
%    Chi_q_PowerLaw_fit_Voro
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
%
%-

chiq_path='E:\Chieco\Hyperuniformity\HU foams\fourier space data\voronoi data';

type={'Poisson'};%;'Einstein_square_delta026'%
pattern_vec={'pointsINI_points_weighted'};%{'centroid_points_weighted'};,

color_mat=[194,0,0;204,136,0;0,0,204]/255;
pl_color=[0,255,0;0,255,255]/255;

for i1=numel(pattern_vec):numel(pattern_vec)
    pattern=pattern_vec{i1};    
    for i2=1:numel(type)        
        chiq=dlmread([chiq_path '\Chi_q ' type{i2} ' ' pattern ' N_cells_ 500000 L_sys_32768.txt']);
        chiq_fit=chiq(2:end,:);
        
        x_all=chiq_fit(:,1);
        if i2==1
            x_min=0.0002;%min(x_all);
            x_max=0.3;%10*x_min;
        else
            x_min=0.001;%min(x_all);
            x_max=0.03;%10*x_min;
        end
            
        y_all=chiq_fit(:,3);
        y_err_all=chiq_fit(:,4);
        
        x_fit_initial=x_all(and(x_all>=x_min,x_all<=x_max));
        y_fit_initial=y_all(and(x_all>=x_min,x_all<=x_max));
        y_err_fit_initial=y_err_all(and(x_all>=x_min,x_all<=x_max));
        
        x_fit=x_fit_initial(y_fit_initial>=6E-10);%x_fit_initial;%
        y_fit=y_fit_initial(y_fit_initial>=6E-10);%y_fit_initial;%
        y_err_fit=y_err_fit_initial(y_fit_initial>=6E-10);
        y_weights=1./(y_err_fit.^2);
        
%         pl_fit=@(coeffs,x_in)coeffs(1)*x_in.^3.63;
        pl_fit=@(coeffs,x_in)coeffs(1)*x_in.^coeffs(2);
         
        coeffs_in=[10,3.5];
        [my_fit,Rs,Jac,CovB,MSE]=nlinfit(x_fit,y_fit,pl_fit,coeffs_in,'Weights',y_weights);        
        coeffs_err=nlparci(my_fit,Rs,Jac);
        my_err=max(abs([my_fit(2)-coeffs_err(2),my_fit(2)-coeffs_err(2)]));
        pow=my_fit(2);%3.7%
        [my_fit,my_err]
        
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
        
%         filewrite_name=[chiq_path '\power law fits\Chi_power_law ' type{i2} ' ' pattern ' N_cells_500000 L_sys_32768 epsioln=' num2str(pow) '.txt'];
%         dlmwrite(filewrite_name,[x_plot',y_plot'],'newline','pc')
    end
end
plotedit on