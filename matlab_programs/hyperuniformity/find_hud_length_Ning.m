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

window='circle';

edit=['variance_avg_s1E5'];
h_edit=['hu length'];
path='D:\Chieco\Hyperuniformity\HU defects\lattice triangular data';
position_path=[path '\particle positions\N_initial=1E6'];
variance_path=[path '\real space averages\variance'];
h_path=[path '\real space averages\hu length'];
h_fit_path=[path '\real space averages\hu length\Einstein swap fits'];


%These are the values of phi we are investigating
 my_num=[0];%,0.00035,0.00065,0.0008,0.001,0.002];%,0.2,0.5];%,1,1.5,2,5,10];

%  my_num_str={'point_5 vnum_point_0','point_5 vnum_point_001'};
 my_num_str={'point_5 swap_NtotXpoint_0','point_5 swap_NtotXpoint_01','point_5 swap_NtotXpoint_02',...
             'point_5 swap_NtotXpoint_03','point_5 swap_NtotXpoint_05','point_5 swap_NtotXpoint_1',...
             'point_5 swap_NtotXpoint_3','point_5 swap_NtotXpoint_5'};%
 p_vec=[0,0.01,0.02,0.03,0.05,0.1,0.3,0.5];
 c_vec=zeros(numel(p_vec),3);
 pattern=my_num_str;
 L_sys=1;
 A_sys=L_sys^2;
 num_runs=1;
 
 i2_start=1 ; i2_fin=numel(my_num_str) ; step=1; count=1;
 for i1=1:numel(pattern)
     for loop=1:num_runs
         lp_text=num2str(loop);
         f_id=pattern{i1};
         xyr_plus=dlmread([position_path '\loop_1 Einstein_strict_swap delta_' my_num_str{i1} '.txt']);
         tag=[window ' Einstein_strict_swap delta_'] ;
         file_name=['\' edit ' ' tag my_num_str{i1} '.txt'];
         path_to_read=[variance_path '\circular window data' file_name];
         variance_data=dlmread(path_to_read);
         particle_area=variance_data(1,3);
         my_phi=sum(pi*(xyr_plus(2:end,3).^2))/A_sys;
         %my_phi=my_num(i2)/10^(ceil(log10(my_num(i2)+1)));%my_phi=pi*sum(xyr(:,3).^2)/A_sys;%my_num(i2,loop)*particle_area/A_sys;%
         %some data is held differently, this is just so I can change them to whatever format
         %I need and not have to change any variables
         lat_sp=variance_data(1,2);
         sqrt_avga=sqrt(variance_data(1,3));
         length_norm=sqrt_avga;%lat_sp;%
         variance_in=[variance_data(:,1)*lat_sp,variance_data(:,6)];
         err_calc=zeros(numel(variance_in(:,1)),4);
         %Now we need to calculate the error for the variance in phi
         %for continuum model we have S=[1-(1-q)^M ]/q   where  q=Vwin/Vsys
         %but for now V_sys=1
         for err=1:numel(err_calc(:,1))
             qs=variance_in(err,1)^2;
             S_num=(1-(1-qs)^variance_data(err,5))/qs;
             err_small=(variance_data(err,4)/(variance_in(err,1)^6))/(S_num*my_phi);
             err_big=(2*variance_in(err,2)^2)/(S_num-1);
             err_calc(err,1)=sqrt(err_small+err_big);
             err_calc(err,2:4)=[qs,S_num,variance_in(err,1)/length_norm];
         end
         variance_in=[variance_in,err_calc(:,1)];
         file_write=[h_path '\' h_edit ' ' tag my_num_str{i1} '.txt'];
         file_fit_write=[h_fit_path '\fit_alpha ' h_edit ' ' tag my_num_str{i1} '.txt'];
         if strcmp(window,'square')~=1
             [hs,error]=calculate_h_circle(variance_in,particle_area);
             loglog(2*hs(:,1)/length_norm,hs(:,2)/(length_norm))                 
             hold on
             dlmwrite(file_write,[2*hs(:,1)/length_norm,hs(:,2)/length_norm,error/length_norm], 'delimiter' ,',', 'newline','pc');
         else
             [hs,error]=calculate_h_square(variance_in,particle_area);
             loglog(hs(:,1)/length_norm,hs(:,2)/length_norm,'LineWidth',2)
             hold on
             ell_s=5.99*length_norm;
             ell_f=L_sys/2;
             h_fit=[hs(and(hs(:,1)>=ell_s,hs(:,1)<=ell_f),:),error(and(hs(:,1)>=ell_s,hs(:,1)<=ell_f))/12];
             if i1==1                 
                 h_0=mean(h_fit(:,2));  
                 h_0_fit=h_fit;
                 h_0_of_L=hs;
                 plot_xs=h_0_of_L(:,1);
                 plot_ys_lin=zeros(numel(hs(:,1)),1)+h_0;
                 plot(plot_xs/sqrt_avga,plot_ys_lin/sqrt_avga,'LineWidth',2)
                 c_vec(i1)=h_0;
             else
                 coeffs_in=[0.001];
                 eff=p_vec(i1);
                 y_weights=1./(h_fit(:,3).^2);
                 h_fit_lin=@(coeffs,ell) h_0_fit(:,2)+coeffs(1)*ell/4;
                 [coeffs_lin,Rs,Jac,CovB,MSE]=nlinfit(h_fit(:,1),h_fit(:,2),h_fit_lin,coeffs_in,'Weights',y_weights);
                 coeffs_err=nlparci(coeffs_lin,Rs,Jac);
                 my_err=2*max(abs([coeffs_lin-coeffs_err(:,1),coeffs_lin-coeffs_err(:,1)]),[],2);
                 % h_fit_sqrt=@(coeffs,ell) (ell/2).*(1-sqrt(1-2*h_0_fit(:,2)./ell).^2-coeffs(1));
                 % coeffs_sqrt=nlinfit(h_fit(:,1),h_fit(:,2),h_fit_sqrt,coeffs_in,'Weights',y_weights);
                 plot_xs=h_0_of_L(:,1);
                 plot_ys_lin=h_0_of_L(:,2)+coeffs_lin(1)*plot_xs/4;
                 %                  plot_ys_sqrt=(plot_xs/2).*(1-sqrt(1-2*h_0_of_L(:,2)./plot_xs).^2-coeffs_sqrt(1));
                 %                  plot(plot_xs/sqrt_avga,plot_ys_sqrt/sqrt_avga,'LineWidth',1.5,'LineStyle', ':');
                 plot(plot_xs/sqrt_avga,plot_ys_lin/sqrt_avga,'LineWidth',1.5,'LineStyle', '--');
                 c_vec(i1,:)=[eff,coeffs_lin/eff,my_err/eff];
                 [p_vec(i1),coeffs_lin*eff]
                 bark=0;
             end
             hold on
             dlmwrite(file_write,[hs(:,1)/length_norm,hs(:,2)/length_norm,error/length_norm], 'delimiter' ,',', 'newline','pc');
             dlmwrite(file_fit_write,[plot_xs/length_norm,plot_ys_lin/length_norm], 'delimiter' ,',', 'newline','pc');
         end
     end
 end
 c_vec=c_vec(2:6,:);
 
% %  coeffs_in=[0,(1/1.4)^2];
% %  
% %  del_fit=@(deltas,eff)(deltas(2)*eff-deltas(1)^2*eff.^2)./(1+deltas(1)*eff);
% %  del_full_fit=@(deltas,eff)eff*(deltas(1)+deltas(2))./(1+deltas(1)*eff);
% %  
% %  coeffs_del=nlinfit(p_vec(2:end-2)',c_vec(2:end-2,3),del_fit,coeffs_in)
% %  coeffs_del_full=nlinfit(p_vec(2:end)',c_vec(2:end,3),del_full_fit,coeffs_in)
% %  hold off
% %  plot(p_vec(2:end)',c_vec(2:end,3)','or')
% %  plot_xs=(0:0.005:1);
% %  plot_ys_del=(coeffs_del(2)*plot_xs-(coeffs_del(1)*plot_xs).^2)./(1+coeffs_del(1)*plot_xs);
% %  hold on
% %  plot(plot_xs,plot_ys_del)
% %  plot_ys_del_full=plot_xs*(coeffs_del_full(1)+coeffs_del_full(2))./(1+coeffs_del_full(1)*plot_xs);
% %  plot(plot_xs,plot_ys_del_full)
%  
% %  delta1=a_big/a_small-1;
% %  delta2=a_small/a_big-1;
% %  avg_delta=(delta1+delta2)/2;
% %  avg_delta_sq=(delta1^2+delta2^2)/2;
% %  [avg_delta,avg_delta_sq]
% %   
% %  for i2=1:numel(p_vec)
% %      eff=p_vec(i2);
% %      plot_xs=h_0_of_L(:,1);
% %      plot_ys_del=(1+coeffs_del(1)*0)*h_0_of_L(:,2)+...
% %          ((coeffs_del(2)*eff-coeffs_del(1)^2*eff.^2)/(1+coeffs_del(1)*eff))*plot_xs/4;
% %      plot_ys_del_calc=(1+avg_delta*0)*h_0_of_L(:,2)+...
% %          ((avg_delta_sq*eff-avg_delta^2*eff.^2)/(1+avg_delta*eff))*plot_xs/4;
% % %      plot_ys_del_full=plot_xs/2.*(1-sqrt((1+coeffs_del_full(1)*0)*(1-2*h_0_of_L(:,2)./plot_xs).^2-...
% % %          eff*(coeffs_del_full(1)+coeffs_del_full(2))/(1+coeffs_del_full(1)*eff))); 
% %      plot(plot_xs/sqrt_avga,plot_ys_del/sqrt_avga,'LineWidth',1.5,'LineStyle', ':')
% %      plot(plot_xs/sqrt_avga,plot_ys_del_calc/sqrt_avga,'LineWidth',1.5,'LineStyle', '-.')
% % %      plot(plot_xs/sqrt_avga,plot_ys_del_full/sqrt_avga,'LineWidth',1.5,'LineStyle', '-.')
% %      file_fit_write=[h_fit_path '\fit_del ' h_edit ' ' tag my_num_str{i2} '.txt'];
% %      dlmwrite(file_fit_write,[plot_xs/length_norm,plot_ys_del/length_norm], 'delimiter' ,',', 'newline','pc');
% %      file_fit_write=[h_fit_path '\fit_del_calc ' h_edit ' ' tag my_num_str{i2} '.txt'];
% %      dlmwrite(file_fit_write,[plot_xs/length_norm,plot_ys_del_calc/length_norm], 'delimiter' ,',', 'newline','pc');
% %  end
 
 
 'all done'
 
  