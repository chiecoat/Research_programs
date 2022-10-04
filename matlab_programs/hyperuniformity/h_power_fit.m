% +
% NAME: h_power_fit
%
% PURPOSE:
%     This program will find find a power law that h(L)abides by. It will
%     go something like h~L^(1-epsilon), where epsilon usually falls
%     between 0<=epsilon<=1. We've found that for some patterns epsilon is
%     greater than 1.
%
% CATEGORY:
%     Hyperuniformity
%
% CALLING SEQUENCE:
%    h_min_he
%
% INPUTS: None at the moment as it is a stand alone program
%
% OPTIONAL INPUTS: (none)
%
% KEYWORD PARAMETERS: (none)
%
% OUTPUTS: Will write a text, that has [L/sqrt(<a>),h_power_fit]

% SIDE EFFECTS: (none)
%
% MODIFICATION HISTORY:
%    written by: A. Chieco, UPenn, September 2016
%
%-
clear all
close all

path='D:\Chieco\Hyperuniformity\HU foams';

h_path=[path '\real space averages\hu length\foam data'];
h_fit_path=[path '\real space averages\hu length\foam data\power law fits'];
h_edit='\hu length circle_windows' ;

%These are the values of phi we are investigating
pattern_type={'arc_points_weighted'};
my_nums={'500','250','125'};
foam_color_mat=[190,190,190;120,120,120;60,60,60;]/255;

i1_start=1 ; i1_fin=numel(pattern_type) ; step=1;

%We want to fit h=K*L^(1-epsilon)
h_power=@(coeffs,xdata)(coeffs(1)*xdata.^(1-coeffs(2)));
coeffs_mat=zeros(numel(my_nums),2);
for i1=i1_start:step:i1_fin
    for i2=1:numel(my_nums)
        file_read=[h_path h_edit ' '  pattern_type{i1} ' imi_' my_nums{i2} '_avg.txt'];
        hs=dlmread(file_read);
        %This is the particle size, the particles are monodisperse
        %         avg_a=1/num_vec(i2); sqrt_avg_a=sqrt(avg_a);
        hs_for_fit=hs(and(hs(:,1)>1,hs(:,1)<=max(hs(:,1))*0.3/0.8),1:3);
        max(hs(:,1))*0.3/0.8
        ell=hs_for_fit(:,1);
        hs_in=hs_for_fit(:,2);
        del_hs=hs_for_fit(:,3);
        y_weights=1./(del_hs.^2);
        coeffs_in=[1,0.5];
        h_fit_coeffs=nlinfit(ell,hs_in,h_power,coeffs_in,'Weights',y_weights)
        coeffs_mat(i2,:)=h_fit_coeffs;
        x_plot=logspace(-1,log10(max(ell)),500);
        h_to_plot=h_fit_coeffs(1)*x_plot.^(1-h_fit_coeffs(2));
        loglog(hs(:,1),hs(:,2),'LineWidth',2,'Color',foam_color_mat(i2,:));
        hold on
        plot(x_plot,h_to_plot,'--r','LineWidth',2);
        file_write=[h_fit_path '\h_power_law '  pattern_type{i1} ' imi_' my_nums{i2} '_avg epsilon=' num2str(h_fit_coeffs(2)) '.txt'];
        dlmwrite(file_write,[x_plot',h_to_plot'],'delimiter',',','newline','pc');
        if i2==1
            data_keep=[ell,hs_in];
        else
            data_ys=interp1(ell,hs_in,data_keep(:,1));
            data_keep=[data_keep(:,1),data_keep(:,2)+data_ys];
        end
    end
end
coeffs_in=[mean(coeffs_mat(:,1)),mean(coeffs_mat(:,2))];
data_keep(:,2)=data_keep(:,2)/3;
plot(data_keep(:,1),data_keep(:,2),'r')
[h_fit_coeffs,Rs,Jac,CovB,MSE]=nlinfit(data_keep(:,1),data_keep(:,2),h_power,coeffs_in)
coeffs_err=nlparci(h_fit_coeffs,Rs,Jac);
% my_err=max(abs([h_fit_coeffs(2)-coeffs_err(2),h_fit_coeffs(2)-coeffs_err(2)]));
% my_err=3*my_err;
my_err=0.15%2*max(abs(h_fit_coeffs(2)-coeffs_mat(:,2)));
my_coeffs_2=h_fit_coeffs(2);
[h_fit_coeffs(2),my_err]

x_plot=logspace(-1,2,500);
h_tot_plot=h_fit_coeffs(1)*x_plot.^(1-h_fit_coeffs(2));
plot(x_plot,h_tot_plot,'--k','LineWidth',2);

y_plot_err_plus=h_fit_coeffs(1)*x_plot.^(1-h_fit_coeffs(2)+my_err);
loglog(x_plot,y_plot_err_plus,'--','Color',[60 60 60]/255,'LineWidth',2);

y_plot_err_minus=h_fit_coeffs(1)*x_plot.^(1-h_fit_coeffs(2)-my_err);
loglog(x_plot,y_plot_err_minus,'--','Color',[60 60 60]/255,'LineWidth',2);

plot(logspace(-4,4,500),logspace(-4,4,500)/2,'--r','Linewidth',2)

file_write=[h_fit_path '\h_power_law '  pattern_type{i1} ' tot_avg epsilon=' num2str(h_fit_coeffs(2)) '.txt'];
dlmwrite(file_write,[x_plot',h_tot_plot'],'delimiter',',','newline','pc');

keyboard

%we are trying to see if the power law fit to the h(L) data is the same
%power for the spectral density data. Here we fit each spectral density
%curve with the same power and with a different power and see which works.
chi_path=[path '\fourier space data\foam data'];
chi_fit_path=[path '\fourier space data\foam data\power law fits'];
chi_edit='\Chi_q_avgellnorm' ;

i1_start=1 ; i1_fin=numel(pattern_type) ; step=1;

%We want to fit h=K*L^(1-epsilon)
chi_power_lock=@(coeffs,xdata)(coeffs(1)*xdata.^(my_coeffs_2));
chi_power_float=@(coeffs,xdata)(coeffs(1)*xdata.^(coeffs(2)));
coeffs_mat_chi=zeros(numel(my_nums),4);
for i1=i1_start:step:i1_fin
    for i2=1:numel(my_nums)
        file_read=[chi_path chi_edit ' '  pattern_type{i1} ' width_eq_1 square_imi_' my_nums{i2} ' L_sys_2820.txt'];
        chi_q=dlmread(file_read);
        %This is the particle size, the particles are monodisperse
%         avg_a=1/num_vec(i2); sqrt_avg_a=sqrt(avg_a);
        chi_q_for_fit=chi_q(and(chi_q(:,1)<0.3,chi_q(:,1)>0),:);
        qs=chi_q_for_fit(:,1);
        chi_q_in=chi_q_for_fit(:,3);
%         del_chi_q=chi_q_for_fit(:,3);
        y_weights=1./(sqrt(chi_q_for_fit(:,3).^2));
        coeffs_lock=[1];
        chi_lock_coeffs=nlinfit(qs,chi_q_in,chi_power_lock,coeffs_lock,'Weights',y_weights);
        coeffs_mat_chi(i2,1:2)=[chi_lock_coeffs,my_coeffs_2];
        x_plot=logspace(log10(min(qs)/10),log10(max(qs)),500);
        chi_to_plot=chi_lock_coeffs(1)*x_plot.^(my_coeffs_2);
        loglog(chi_q(:,1),chi_q(:,3),'LineWidth',2);
        hold on
        plot(x_plot,chi_to_plot,'--r','LineWidth',2);
        file_write=[chi_fit_path '\chi_power_law '  pattern_type{i1} ' square_imi_' my_nums{i2} ' L_sys_2820 epsilon=' num2str(my_coeffs_2) '.txt'];
        dlmwrite(file_write,[x_plot',chi_to_plot'],'delimiter',',','newline','pc');
        coeffs_float=[1,my_coeffs_2];
        chi_float_coeffs=nlinfit(qs,chi_q_in,chi_power_float,coeffs_in,'Weights',y_weights);
        coeffs_mat_chi(i2,3:4)=chi_float_coeffs;
        x_plot=logspace(log10(min(qs)/10),log10(max(qs)),500);
        chi_to_plot=chi_float_coeffs(1)*x_plot.^(chi_float_coeffs(2));
        plot(x_plot,chi_to_plot,'--b','LineWidth',2);
        if i2==1
            data_keep=[qs,chi_q_in];
        else
            data_ys=interp1(qs,chi_q_in,data_keep(:,1));
            data_keep=[data_keep(:,1),data_keep(:,2)+data_ys];
        end
    end
end

my_chi_coeffs_1=mean(coeffs_mat_chi(:,1));
chi_fixed_power=nlinfit(qs,chi_q_in,chi_power_lock,my_chi_coeffs_1,'Weights',y_weights);
x_plot=logspace(-4,2,500);
chi_tot_plot=chi_fixed_power*x_plot.^(my_coeffs_2);
plot(x_plot,chi_tot_plot,'--k','LineWidth',2);

file_write=[chi_fit_path '\chi_power_law '  pattern_type{i1} ' tot_avg epsilon=' num2str(my_coeffs_2) '.txt'];
dlmwrite(file_write,[x_plot',chi_tot_plot'],'delimiter',',','newline','pc');
% 
'all done'