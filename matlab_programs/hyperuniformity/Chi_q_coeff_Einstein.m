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

chiq_path='F:\Chieco\Hyperuniformity\HU defects\area defects\lattice square data\fourier space data';
chiq=dlmread([chiq_path '\centralpoint Chi_q_Einstein_delta_point_0.31188 f=0 ad_mag=0 phi_1_1 N_sys_1E7 L_sys_8192.txt']);
chiq=chiq(2:end-1,:);

chiq_a=chiq(chiq(:,1)>1.1,:);

beta1=1;

alpha1_avg_vec=2*chiq_a(:,2)/(beta1*pi);
alpha1_avg=mean(alpha1_avg_vec);

ell1=beta1*pi./chiq(:,1);
h1_from_chiq=chiq(:,2)./(alpha1_avg*chiq(:,1)); 

loglog(ell1,h1_from_chiq);
hold on

beta2=2;

alpha2_avg_vec=2*chiq_a(:,2)/(beta2*pi);
alpha2_avg=mean(alpha2_avg_vec);

ell2=beta2*pi./chiq(:,1);
h2_from_chiq=chiq(:,2)./(alpha2_avg*chiq(:,1)); 

loglog(ell2,h2_from_chiq);
 
h_path='F:\Chieco\Hyperuniformity\HU low discrepancy patterns\real space averages\hu length';
h_data=dlmread([h_path '\hu_length Einstein_square_N1E7 phi_1_avg delta_point_0.31188f=0 ad_mag=0 calc_from_circ_short.txt']);


loglog(h_data(:,1),h_data(:,2))
loglog(h_data(:,1),h_data(:,1)/2,'Linewidth',2,'Color','r')
loglog([1,1E4],[mean(h_data(h_data(:,1)>10,2)),mean(h_data(h_data(:,1)>10,2))],'-.m','Linewidth',2)    

bark=0;





