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
clear all
close all

edit='variance Einstein_square_N1E7';
h_edit='hu_length Einstein_square_N1E7';

path='F:\Chieco\Hyperuniformity\HU low discrepancy patterns\real space averages';
variance_path=[path '\variance'];
h_path=[path '\hu length'];

%These are the values of nosie we are investigating
eff_val=[0];
eff_titles={'0'};
kick_vec=[0.31188];
kick_titles={'0.31188'};

i2_start=1 ; i2_fin=numel(eff_titles) ; step=1;
 for i2=i2_start:step:i2_fin
   eff=eff_titles{i2};
   kick=kick_titles{i2};
   my_phi=1;
   txt_phi=num2str(my_phi);
   file_name=['\' edit ' ' 'phi_' txt_phi '_avg delta_point_' kick ' f=' eff ' ad_mag=0 calc_from_circ_short.txt'];
   path_to_read=[variance_path file_name];
   variance_data=dlmread(path_to_read);
   particle_area=variance_data(1,2);
   p=find(variance_data(:,5) > 0);
   variance_data=[variance_data(p,1),variance_data(p,5),variance_data(p,5)];
   keep=find(variance_data(:,2) > 0);
  %some data is held didfferenty, this is just so I can change them to whatever format
  %I need and not have to change any variables
   sqrt_avga=sqrt(particle_area(1));
   variance_in=[variance_data(:,1)*sqrt_avga,variance_data(:,2:3)];
   [hs,error]=calculate_h_circle(variance_in,particle_area(1));
   file_write=['\' h_edit ' ' 'phi_' txt_phi '_avg delta_point_' kick 'f=' eff ' ad_mag=0 calc_from_circ_short.txt'];
   dlmwrite([h_path file_write],[2*hs(:,1)/sqrt_avga,hs(:,2)/sqrt_avga,error/sqrt_avga],'newline','pc');
   bark=0;
 end
'all done'

