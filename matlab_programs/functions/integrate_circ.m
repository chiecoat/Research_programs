% +
% NAME: foam_track_2d
%
% PURPOSE: We numerically integrate coarsening data based on circularity
%
% CATEGORY:
%     Foam coarsening
%
% CALLING SEQUENCE:
%    integrate_circ(time,[areas,areas_err],[circ,circ_err],n_s,r_curv,k_prime)
%
% INPUTS: 
%    
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
%    written by: A. Chieco, UPenn, Feb. 2018
function areas_expected=integrate_circ(time,area_info,circ_info,elong_info,A_start,n_s,r_curv,k_prime,gap)

n_frames=numel(time);
VN_ex=n_s-6;
A_0=A_start;

areas_expected=zeros(n_frames,1);
area_smooth=smoothdata(area_info(:,1),'movmean',numel(area_info(:,1))/10);
circ_smooth=smoothdata(circ_info(:,1),'movmean',numel(circ_info(:,1))/10);
elong_smooth=smoothdata(elong_info(:,1),'movmean',numel(elong_info(:,1))/10);

for i1=1:n_frames
    if i1==1
        areas_expected(i1)=A_0;
    else
        beta=6*n_s*r_curv/sqrt(3*pi*area_smooth(i1-1,1));
        denom=1-(1-pi/4)*(sqrt(4*pi)*elong_smooth(i1-1,1)*(r_curv^2)/(gap*sqrt(area_smooth(i1-1,1))));
        c_int=trapz(time(i1-1:i1),circ_smooth(i1-1:i1));
        areas_expected(i1)=areas_expected(i1-1)+(k_prime/denom)*(VN_ex*(time(i1)-time(i1-1))+beta*c_int);        
    end   
end




end







%-