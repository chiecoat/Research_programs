% +
% NAME: find_h_length
%
% PURPOSE:
%     This program will solve for h, the "hyperuniformity length" from the variance
%     of a point pattern.
%
% CATEGORY:
%     Hyperuniformity
%
% CALLING SEQUENCE:
%    find_h_length
%
% INPUTS: None at the moment as it is a stand alone program
%
% OPTIONAL INPUTS: (none)
%
% KEYWORD PARAMETERS: (none)
%
% OUTPUTS: Will write a text, that has h computed for different point patterns
%          read in

% SIDE EFFECTS: (none)
%
% MODIFICATION HISTORY:
%    written by: A. Chieco, UPenn, January  2015
%
%-
% clear all
% close all

path='F:\Chieco\Hyperuniformity\HU Jammed Packings from Ning Xu\nearly phi_c\tol_tests';
variance_path=[path '\loop data\variance'];
variance_avg_path=[path '\real space averages\variance\sum averages'];
vrat_avg_path=[path '\real space averages\variance ratio\sum averages'];

pattern='square_N100000_tol1.0e-8_s1E5';

var_edit=['\variance ' pattern];

vrat_edit=['\var_ratio ' pattern];

h_path=[path '\loop data\hu length'];
h_path_gen=[path '\real space averages\hu length\sum averages'];
h_avg_path=h_path_gen;
h_edit=['\hu length ' pattern ];

%These are the values of phi we are investigating
my_phis=[8412];
info_vec=zeros(numel(my_phis),2);

i2_start=1 ; i2_fin=numel(my_phis) ; step=1;
n_runs=5;
n_loops_vec=zeros(1,numel(my_phis))+n_runs;
% n_configs=n_loops*i2_fin;
n_configs=sum(n_loops_vec);
above_jamming_data=[0];

for i2=i2_start:step:i2_fin
    my_phi=my_phis(i2);%/(10^ceil(log10(my_phis(i2)+1)));
    txt_phi=num2str(my_phis(i2));
    n_loops=n_loops_vec(i2);
    for loop=1:n_loops
        lp_text=num2str(loop);
        file_name=[var_edit ' ' txt_phi '_loop_' lp_text '.txt'];
        path_to_read=[variance_path file_name];
        variance_data=dlmread(path_to_read);
        particle_area=variance_data(1,2);
        %some data is held didfferenty, this is just so I can change them to whatever format
        %I need and not have to change any variables
        sqrt_avga=sqrt(particle_area(1));
        variance_in=[variance_data(:,1)*sqrt_avga,variance_data(:,5),zeros(numel(variance_data(:,1)),1)];
        sqrt_avga=sqrt(particle_area);
        err_calc=zeros(numel(variance_in(:,1)),4);
        %Now we need to calculate the error for the variance in phi
        %for continuum model we have S=[1-(1-q)^M ]/q   where  q=Vwin/Vsys
        for err=1:numel(err_calc(:,1))
            qs=variance_in(err,1)^2;
            S_num=(1-(1-qs)^variance_data(err,4))/qs;
            err_small=(variance_data(err,3)/variance_in(err,1)^6)/(S_num*my_phi);
            err_big=(2*variance_in(err,2)^2)/(S_num-1);
            err_calc(err,1)=sqrt(err_small+err_big);
            err_calc(err,2:4)=[qs,S_num,variance_in(err,1)/sqrt_avga];
        end
        variance_in(:,3)=err_calc(:,1);
        [hs,error]=calculate_h_square(variance_in,particle_area(1));
        file_write=[h_path h_edit ' ' txt_phi '_loop_' lp_text '.txt'];
        dlmwrite(file_write,[hs(:,1)/sqrt_avga,hs(:,2)/sqrt_avga,error/sqrt_avga],'Newline','pc');
        if loop==1
            variance_sum=variance_in;
        else
            variance_sum=[variance_in(:,1),variance_sum(:,2)+variance_in(:,2),...
                variance_sum(:,3)+variance_in(:,3)];
        end
        %         if and(i2==i2_start,loop==1)==1
        %             h_above_jamming_data=[hs(:,1),hs(:,2)/sqrt_avga,error/sqrt_avga];
        %         else
        %             h_above_jamming_data=[hs(:,1),h_above_jamming_data(:,2)+hs(:,2)/sqrt_avga,...
        %                 h_above_jamming_data(:,3)+error/sqrt_avga];
        %         end
    end
    variance_avg=[variance_sum(:,1),variance_sum(:,2)/n_loops,variance_sum(:,3)/(sqrt(n_loops)*n_loops)];
    file_var_avg=[variance_avg_path var_edit ' ' txt_phi '_avg.txt'];
    dlmwrite(file_var_avg,[variance_avg(:,1)/sqrt_avga,variance_avg(:,2:3)],'Newline','pc');
    vrat_avg=[variance_sum(:,1),(variance_sum(:,1).^2).*variance_sum(:,2)/(n_loops*particle_area),(variance_sum(:,1).^2).*variance_sum(:,3)/(n_loops*particle_area)];
    file_vrat_avg=[vrat_avg_path vrat_edit ' ' txt_phi '_avg.txt'];
    dlmwrite(file_vrat_avg,[vrat_avg(:,1)/sqrt_avga,vrat_avg(:,2:3)],'Newline','pc');
    %     %here we collect all of the data above jamming into one matrix
    %     if and(my_phi>0.84,numel(above_jamming_data==1))==1
    %         above_jamming_data=[variance_avg(:,1),variance_avg(:,2)*my_phi,...
    %             variance_avg(:,3)*my_phi];
    %     end
    %     if and(my_phi>0.84,numel(above_jamming_data>1))==1
    %         above_jamming_data=[variance_avg(:,1),above_jamming_data(:,2)+variance_avg(:,2)*my_phi,...
    %             above_jamming_data(:,3)+variance_avg(:,3)*my_phi];
    %     end
    info_vec(i2,:)=[my_phi,particle_area(1)];
    [h_avg,error_calc]=calculate_h_square([variance_avg(:,1),variance_avg(:,2),variance_avg(:,3)],particle_area(1));
    error_avg=error_calc;
    err_elim1=3*(1:numel(error_avg)/3);
    err_elim2=err_elim1-1;
    error_avg([err_elim1(1:end-1)';err_elim2'])=0;
    file_write=[h_avg_path h_edit ' ' txt_phi '_avg.txt'];
    dlmwrite(file_write,[h_avg(:,1)/sqrt_avga,h_avg(:,2)/sqrt_avga,error_avg/(sqrt(n_loops)*sqrt_avga)],'Newline','pc');
    bark=0;
    %     %This fits the initial linear deviation%%%%%%%%%%%%%%%%%%%
    %     h_fit=[h_avg(:,1)/sqrt_avga,h_avg(:,2)/sqrt_avga,error_calc/(12*sqrt_avga)];
    %     ell_min=10; ell_max=L_c_vec(i2);
    %     h_keep=find(and(h_fit(:,1)>=ell_min,h_fit(:,1)<=ell_max));
    %     ell=h_fit(h_keep,1);
    %     hs=h_fit(h_keep,2);
    %     del_hs=h_fit(h_keep,3);
    %     y_weights=1./(del_hs.^2);
    %     %we fit to a function h(L)=h_e+beta*L/4
    %     h_e=0.084;
    %     coeffs_in=[0.002];
    %     hfit_forLe=@(coeffs,l_fit)h_e+coeffs(1)*l_fit/4;
    %     [coeffs_lin,Rs,Jac,CovB,MSE]=nlinfit(ell,hs,hfit_forLe,coeffs_in,'Weights',y_weights);
    %     coeffs_err=nlparci(coeffs_lin,Rs,Jac);
    %     my_beta_err=max(abs([coeffs_lin(1)-coeffs_err(1,:),coeffs_lin(1)-coeffs_err(1,:)]),[],2);
    % %     my_h_err=2*max(abs([coeffs_lin(2)-coeffs_err(2,:),coeffs_lin(2)-coeffs_err(2,:)]),[],2);
    %     L_fin=1E4;
    %     plot_xs=(1.1:(L_fin-1.1)/1E3:L_fin);
    %     my_h_fit=h_e+coeffs_lin(1)*plot_xs/4;
    %     loglog(h_fit(:,1),h_fit(:,2),'LineWidth',1.5)
    %     hold on
    %     loglog(plot_xs,my_h_fit,'--','LineWidth',1.5)
    %     [coeffs_lin(1),my_beta_err]
    %     h_file_fit=['\' pattern ' hu length fit_he_beta'];
    %     file_write=[h_fit_path h_file_fit ' phi_' txt_phi '.txt'];
    %     dlmwrite(file_write,[plot_xs',my_h_fit'],'Newline','pc');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     loglog(h_keep(:,1)/sqrt_avga,h_swath(:,2)/sqrt_avga,'b','LineWidth',1)
    %     hold on
    %     loglog(h_keep(:,1)/sqrt_avga,h_swath(:,1)/sqrt_avga,'r','LineWidth',1)
    %     loglog(h_keep(:,1)/sqrt_avga,h_swath(:,3)/sqrt_avga,'g','LineWidth',1)
    %     bark=0;
    %     if i2==i2_start
    %         hs_keep=h_avg(:,2)/sqrt_avga;
    %     else
    %         hs_keep=[hs_keep,h_avg(:,2)/sqrt_avga];
    %     end
end
% % avg_a_bar=max(info_vec(:,2));
% % sqrt_avg_a_bar=sqrt(avg_a_bar);
% % above_jamming_data=[above_jamming_data(:,1),...
% %     above_jamming_data(:,2)/(i2_fin*mean(info_vec(:,1))),...
% %     above_jamming_data(:,3)/(i2_fin*mean(info_vec(:,1)))];
% % file_var_avg=[variance_avg_path w_var_edit ' all above_jamming_avg.txt'];
% % dlmwrite(file_var_avg,[above_jamming_data/sqrt(avg_a_bar),above_jamming_data(:,2:3)],'Newline','pc');
% % file_write=[h_avg_path h_edit ' all above_jamming_avg.txt'];
% % %we add a line for a box equal to the system size. The variance is zero, so
% % %h=0 and error=0 as well
% % h_write=find(h_above_jamming_data(:,1)<1E4*sqrt_avg_a_bar);
% % err_elim1=3*(1:numel(error_avg)/3);
% % err_elim2=err_elim1-1;
% % h_above_jamming_data([err_elim1(1:end-1)';err_elim2'],3)=0;
% % h_avg_write=[h_above_jamming_data(h_write,1)/sqrt_avg_a_bar,h_above_jamming_data(h_write,2)/n_configs,...
% %     h_above_jamming_data(h_write,3)/(n_configs*sqrt(n_configs))];
% % dlmwrite(file_write,h_avg_write,'Delimiter',',','Newline','pc');
% 
% h_avg_fit=[h_avg(:,1)/sqrt_avga,h_avg(:,2)/sqrt_avga,error_calc/sqrt_avga];
% % %
% % % % h_DJD_file=[h_above_jamming_data(h_write,1)/sqrt_avg_a_bar,hs_keep,h_above_jamming_data(h_write,2)/n_configs,...
% % % %              h_above_jamming_data(h_write,3)/(n_configs*sqrt(n_configs)),zeros(numel(h_above_jamming_data(h_write,1)),1)];
% % % % for ii=1:numel(h_DJD_file(:,1))
% % % %     err_stdv=std(h_DJD_file(ii,2:4));
% % % %     h_DJD_file(ii,end)=err_stdv;
% % % % end
% % % % file_write=[h_avg_path h_edit ' all_h above_jamming_avg.txt'];
% % % % dlmwrite(file_write,h_DJD_file,'Delimiter',',','Newline','pc');
% % % % bark=0;
% % % %
% % % % % h_avg_fit=[h_above_jamming_data(:,1)/sqrt_avga,h_above_jamming_data(:,2)/n_configs,...
% % % % %     h_DJD_file(:,end)];
% % % % h_e=2*h_e;
% 
% % % % %Now we fit to fine L_e for the data
% L_sys=1;
% ell_min=10; ell_max=300;
% h_e=0.084;
% xs_accept=find(and(h_avg_fit(:,1)>=ell_min,h_avg_fit(:,1)<=ell_max));
% h_keep=find(h_avg_fit(xs_accept,2)<=2*h_e);
% ell=h_avg_fit(xs_accept(h_keep),1);
% h_fit=h_avg_fit(xs_accept(h_keep),2);
% del_hs=h_avg_fit(xs_accept(h_keep),3);
% y_weights=1./(del_hs.^2);
% L_e=300;
% 
% %we fit to a function whose parameters are f, fraction of defects, and delta, average size of defects.
% %We fit the linear regime first by ensuring our y values do not go above a
% %certain percentage of h_e
% coeffs_in=[0.0818,L_e];
% % hfit_forLe=@(coeffs,l_fit)coeffs(1)*(1+l_fit/coeffs(2));
% % coeffs_lin=nlinfit(ell,h_fit,hfit_forLe,coeffs_in,'Weights',y_weights);
% % %We fit again this time over the whole range. What is different is now we
% % %are only fititng for the quadratic coefficient, and we use the values for
% % %h_e and L_e as from the previous fits as granted.
% % coeffs_in=[coeffs_lin(1),coeffs_lin(2),1/coeffs_lin(2)];
% % ell_offset=0;
% % hfit_forL2=@(coeffs,l_fit)coeffs(1)*(1+(l_fit/coeffs(2)).*(1+l_fit*coeffs(3)));
% % % model_quad=fitnlm(ell,h_fit,hfit_forL2,coeffs_in,'Weights',y_weights)
% % % coeffs_quad=model_quad.Coefficients.Estimate;
% % [coeffs_quad,Rs,Jac,CovB,MSE]=nlinfit(ell,h_fit,hfit_forL2,coeffs_in,'Weights',y_weights);
% % coeffs_err=nlparci(coeffs_quad,Rs,Jac);
% % my_err=2*max(abs([coeffs_quad'-coeffs_err(:,1),coeffs_quad'-coeffs_err(:,1)]),[],2);
% %
% % ell_plot_fit=logspace(log10(ell_min),log10(max(ell)));
% % ell_plot_big=logspace(log10(5),log10(5000));
% %
% %
% % del_Le=my_err(2);
% % h_from_fit=coeffs_quad(1)*(1+(ell_plot_fit/coeffs_quad(2)).*(1+ell_plot_fit*coeffs_quad(3)));
% % h_minus=coeffs_quad(1)*(1+(ell_plot_fit/(coeffs_quad(2)-del_Le)).*(1+ell_plot_fit*coeffs_quad(3)));
% % h_plus=coeffs_quad(1)*(1+(ell_plot_fit/(coeffs_quad(2)+del_Le)).*(1+ell_plot_fit*coeffs_quad(3)));
% % h_from_all=coeffs_quad(1)*(1+(ell_plot_big/coeffs_quad(2)).*(1+ell_plot_big*coeffs_quad(3)));
% %
% % loglog(h_avg_fit(:,1),h_avg_fit(:,2),'LineWidth',2,'Color','k')
% % hold on
% % errorbar(h_avg_fit(:,1),h_avg_fit(:,2),h_avg_fit(:,3))
% % plot(ell_plot_big,h_from_all,'LineWidth',2,'Color','m')
% % plot(ell_plot_fit,h_from_fit,'LineWidth',2,'Color','r')
% % plot(ell_plot_fit,h_plus,'LineWidth',2,'Color','g')
% % plot(ell_plot_fit,h_minus,'LineWidth',2,'Color','b')
% %
% % hfit_Le=['\fix_he=' num2str(h_e) ' L_e=' num2str(coeffs_lin(1)) '.txt'];
% % fit_write_f=[h_path_gen '\fits for hu length' hfit_Le];
% % dlmwrite(fit_write_f,[ell_plot'/sqrt_avg_a_bar,h_from_lin'/sqrt_avg_a_bar],'Delimiter',',','newline','pc');
% 
% 
% % he_print=coeffs_quad(1);
% % Lc_print=coeffs_quad(2);
% % k_L=(coeffs_quad(3));
% % beta=he_print/Lc_print;
% % eff=4*(beta)
% bark=0;
% 
% 
% %Now we constrain h(L) to fit to asympotic linear form we calculate from
% %L_e & h_e
% ell=h_avg_fit(and(h_avg_fit(:,1)>=ell_min,h_avg_fit(:,1)<=ell_max),1);
% h_fit=h_avg_fit(and(h_avg_fit(:,1)>=ell_min,h_avg_fit(:,1)<=ell_max),2);
% del_hs=h_avg_fit(and(h_avg_fit(:,1)>=ell_min,h_avg_fit(:,1)<=ell_max),3);
% my_asymptote=0;
% if my_asymptote==0
%     my_new_hs=[ell,h_fit,del_hs];
%     file_id='hu length';
% else
%     my_new_hs=[ell,h_fit,del_hs;max(ell_plot),max(h_from_lin),1E-8];
%     file_id='hu length asymptote';
% end
% new_weights=1./(my_new_hs(:,3).^2);
% 
% 
% %Now we fit only for f, the fraction of vacancies. We will do this for both
% %vacancies, delta=-1. For this we use hfit_softdisk_f
% hfit_softdisk_f=@(coeffs,xdata)(h_e*sqrt(1-coeffs(1))+(xdata/2)*(1-sqrt(1-coeffs(1))));
% %and "added particles", delta=+1. For this we use hfit_softdisk_plus
% hfit_softdisk_plus=@(coeffs,xdata)(xdata/2-(xdata/2).*((1+coeffs(1))*(1-2*h_e./xdata).^2-...
%     2*coeffs(1)/(1+coeffs(1))).^(1/2));
% coeffs_f_in=[0.02];
% %Fit for vacancies
% coeffs_vac=nlinfit(my_new_hs(:,1),my_new_hs(:,2),hfit_softdisk_f,coeffs_f_in,'Weights',new_weights)
% %Fit for added particles
% % coeffs_d_occ=nlinfit(my_new_hs(:,1),my_new_hs(:,2),hfit_softdisk_plus,coeffs_f_in,'Weights',new_weights)
% %plot our findings
% ell_plot=logspace(log10(1),log10(10*ell_max));
% h_from_f=h_e*sqrt(1-coeffs_vac(1))+(ell_plot/2)*(1-sqrt(1-coeffs_vac(1)));
% h_from_f_eye=h_e*sqrt(1-0.02)+(ell_plot/2)*(1-sqrt(1-0.02));
% % h_from_del_plus=(ell_plot/2-(ell_plot/2).*((1+coeffs_d_occ(1))*(1-2*h_e./ell_plot).^2-...
% %                                          2*coeffs_d_occ(1)/(1+coeffs_d_occ(1))).^(1/2));
% loglog(h_avg_fit(:,1),h_avg_fit(:,2))
% hold on
% loglog(ell_plot,h_from_f,'LineWidth',1,'Color','r');
% loglog(ell_plot,h_from_f_eye,'LineWidth',1,'Color','b')
% % plot(ell_plot/sqrt_avg_a_bar,h_from_del_plus/sqrt_avg_a_bar,'LineWidth',1,'Color','m')
% %We writ eout a file here for when delta=-1
% hfit_f=['\vac_input ' file_id ' he=' num2str(h_e/sqrt_avga) ' fit_f=' num2str(coeffs_vac(1)) ' delta=-1.txt'];
% hput_f=['\vac_input ' file_id ' he=' num2str(h_e/sqrt_avga) ' my_f=0.0051 delta=-1.txt'];
% fit_write_f=[h_path_gen '\fits for hu length' hfit_f];
% put_write_f=[h_path_gen '\fits for hu length' hput_f];
% dlmwrite(fit_write_f,[ell_plot',h_from_f'],'Delimiter',',','newline','pc');
% dlmwrite(put_write_f,[ell_plot',h_from_f_eye'],'Delimiter',',','newline','pc');
% %Now we write out the file for when delta=+1
% % hfit_del_plus=['\dubocc_input ' file_id ' he=' num2str(h_e/sqrt_avg_a_bar) ' fit_f=' num2str(coeffs_d_occ(1)) ' delta=1.txt'];
% % fit_write_del_plus=[h_path_gen '\fits for hu length' hfit_del_plus];
% % dlmwrite(fit_write_del_plus,[ell_plot'/sqrt_avga,h_from_f'/sqrt_avga],'Delimiter',',','newline','pc');
% % % %
% % % %
% % % % %Now we fit for both the f and delta data
% % % % hfit_softdisk_vacancies=@(coeffs,xdata)(xdata/2-(xdata/2).*((1+coeffs(1)*coeffs(2))*(1-2*h_e./xdata).^2- ...
% % % %                                          coeffs(1)*coeffs(2)*(1+coeffs(2))/(1+coeffs(1)*coeffs(2))).^(1/2));
% % % %
% % % % delta_minus=-0.15;
% % % % coeffs_df_in=[eff,delta_minus];
% % % % opts=statset('MaxIter',5E4);
% % % % coeffs_vfit=nlinfit(my_new_hs(:,1),my_new_hs(:,2),hfit_softdisk_vacancies,coeffs_df_in,opts,'Weights',new_weights)
% % % % h_from_vfit=ell_plot/2-(ell_plot/2).*((1+coeffs_vfit(1)*coeffs_vfit(2))*(1-2*h_e./ell_plot).^2-...
% % % %                         coeffs_vfit(1)*coeffs_vfit(2)*(1+coeffs_vfit(2))/(1+coeffs_vfit(1)*coeffs_vfit(2))).^(1/2);
% % % % plot(ell_plot/sqrt_avg_a_bar,h_from_vfit/sqrt_avg_a_bar,'LineWidth',1,'Color','b');
% % % % hfit_df=['\vac_fit ' file_id ' vac_he=' num2str(h_e/sqrt_avg_a_bar) ' fit_f=' num2str(coeffs_vfit(1)) ' fit_delta=' num2str(coeffs_vfit(2)) ' .txt'];
% % % % fit_write_df=[h_path_gen '\fits for hu length' hfit_df];
% % % % % dlmwrite(fit_write_df,[ell_plot'/sqrt_avga,h_from_vfit'/sqrt_avga]);
% % % %
% % % % %I want to do this again but with an postive delta input and the same initial f
% % % % %The function may not be well conditioned if we choose delta=+1
% % % % %or delta=-1 as initial guesses so we can rule out double occupancies or
% % % % %vacancies as the defects in the system.
% % % % % hfit_softdisk_doubleocc=@(coeffs,xdata)(xdata/2-(xdata/2).*((1+coeffs(1))*(1-2*h_e./xdata).^2-...
% % % % %                                                 coeffs(1)*(1+coeffs(2))/(1+coeffs(1))).^(1/2));
% % % % % delta_plus=1;
% % % % % coeffs_df_in=[delta_plus*eff,delta_plus];
% % % % % opts=statset('MaxIter',5E4);
% % % % % coeffs_dofit=nlinfit(my_new_hs(:,1),my_new_hs(:,2),hfit_softdisk_doubleocc,coeffs_df_in,opts,'Weights',new_weights)
% % % % % h_from_dofit=ell_plot/2-(ell_plot/2).*((1+coeffs_dofit(1))*(1-2*h_e./ell_plot).^2-...
% % % % %                          coeffs_dofit(1)*(1+coeffs_dofit(2))/(1+coeffs_dofit(1))).^(1/2);
% % % % % plot(ell_plot/sqrt_avg_a_bar,h_from_dofit/sqrt_avg_a_bar,'LineWidth',1,'Color','g');
% % % % % hfit_df=['\dub_occ ' file_id ' he=' num2str(h_e/sqrt_avg_a_bar) ' fit_f=' num2str(coeffs_dofit(1)) ' fit_delta=' num2str(coeffs_dofit(2)) ' .txt'];
% % % % % fit_write_df=[h_path_gen '\fits for hu length' hfit_df];
% % % % % dlmwrite(fit_write_df,[ell_plot'/sqrt_avga,h_from_dofit'/sqrt_avga]);
% % % %
% % % % %One last fit combines the fact that vacancies and additions act in pairs.
% % % % %It assumes the average "size" of a vacancy is zero, but the average value
% % % % %squared is not zero. So we try that fit here. It is equations (15) is
% % % %
% % % % %introduction plus.
% % % % % hfit_softdisk_copair=@(coeffs,xdata)(xdata/2-(xdata/2).*((1+coeffs(1))*(1-2*h_e./xdata).^2 ...
% % % % %                                    -(coeffs(1)+coeffs(2))/(1+coeffs(1))).^(1/2));
% % % % hfit_softdisk_copair=@(coeffs,xdata)(xdata/2).*(1-sqrt(1-4*h_e*(1+coeffs(1)*coeffs(2))*(xdata-h_e)./(xdata.^2)...
% % % %                                     +((coeffs(1)^2)*(coeffs(2)^2)-coeffs(3)*coeffs(1))/(1+coeffs(1)*coeffs(2))));
% % % %
% % % % coeffs_co_in=[coeffs_vfit(1),coeffs_vfit(2),0.25];
% % % % opts=statset('MaxIter',5E4);
% % % % coeffs_cop=nlinfit(my_new_hs(:,1),my_new_hs(:,2),hfit_softdisk_copair,coeffs_co_in,opts,'Weights',new_weights)
% % % % % h_from_copair=(ell_plot/2-(ell_plot/2).*((1+coeffs_copair(1))*(1-2*h_e./ell_plot).^2 ...
% % % % %                                         -(coeffs_copair(1)+coeffs_copair(2))/(1+coeffs_copair(1))).^(1/2));
% % % % h_from_copair=(ell_plot/2).*(1-sqrt(1-4*h_e*(1+coeffs_cop(1)*coeffs_cop(2))*(ell_plot-h_e)./(ell_plot.^2)...
% % % %                                     +((coeffs_cop(1)^2)*(coeffs_cop(2)^2)-coeffs_cop(3)*coeffs_cop(1))/(1+coeffs_cop(1)*coeffs_cop(2))));
% % % % plot(ell_plot/sqrt_avg_a_bar,h_from_copair/sqrt_avg_a_bar,'LineWidth',1,'Color',[1 85/255 0]);
% % % % bark=0;
% % % % % hfit_df=['\copair ' file_id ' vac_he=' num2str(coeffs_vac(2)/sqrt_avg_a_bar) ' fit_f=' num2str(coeffs_copair(1)) ' fit_delta_sq=' num2str(coeffs_copair(2)) ' .txt'];
% % % % % fit_write_df=[h_path_gen '\fits for hu length' hfit_df];
% % % % % dlmwrite(fit_write_df,[ell_plot'/sqrt_avga,h_from_copair'/sqrt_avga]);
% % % %
% % % % % I needed this for Doug to try his own fits
% % % % % h_for_fit_file=[h_path_gen '\hu length for fit.txt'];
% % % % % dlmwrite(h_for_fit_file,my_new_hs,'Delimiter',',','Newline','pc')

