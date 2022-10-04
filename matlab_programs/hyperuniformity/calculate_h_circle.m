% +
% NAME: calculate_h
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
% INPUTS: var_dat,a_p
%
%         var_dat: the first column is window side lengths and the second column is
%                  the variance data from which h is calculate
%
%         a_p: The area of a particle that makes up our paticle packings
%
% OPTIONAL INPUTS: (none)
%
% KEYWORD PARAMETERS: error
%         error: If error is input we calculate +/- values for h
%
% OUTPUTS: lil_h_vec_write
%         lil_h_vec_write: This is a 2xN matrix of window side lengths and hu lengths

% SIDE EFFECTS: (none)
%
% MODIFICATION HISTORY:
%    written by: A. Chieco, UPenn, January  2015
%
%-
function [lil_h_vec_write,error]=calculate_h_circle(var_dat,a_p)

 rad=var_dat(:,1);
%now we calculate h. Because of the nature of small window statistic there may be some L
%values where h might b imaginary. val will help compensate for that
 val=1-(pi*rad.^2/a_p).*var_dat(:,2);
 re=find(val > 0) ; im=find(val < 0);
 im_count=numel(im);
 val1=val(re);
 lil_h_re=rad(re)-sqrt(val1).*rad(re);
 if im_count > 0 
   val2=val(im);
   lil_h_im=rad(im)+sqrt(abs(val2)).*rad(im);
   lil_h_vec=[rad(re),lil_h_re;...
              rad(im),lil_h_im];
   [ell_sort,ell_index]=sort(lil_h_vec(:,1));
   lil_h_vec_write=lil_h_vec(ell_index,:);
 else
   lil_h_vec_write=[rad(re),lil_h_re];
 end

 %here we calculate the error in our hu length, if requested
 if nargout == 2
     %and now we need to repeat this for finding the error in h
     val_max=1-(pi*rad.^2/a_p).*(var_dat(:,2)+var_dat(:,3));
     re_max=find(val_max > 0) ; im_max=find(val_max < 0);
     im_max_count=numel(im_max);
     val1_max=val_max(re_max);
     lil_h_max_re=rad(re_max)-sqrt(val1_max).*rad(re_max);
     if im_max_count > 0 
         val2_max=val_max(im_max);
         lil_h_max_im=rad(im_max)+sqrt(abs(val2_max)).*rad(im_max);
         h_max_vec=[rad(re_max),lil_h_max_re;...
                    rad(im_max),lil_h_max_im];
         [hmax_sort,hmax_index]=sort(h_max_vec(:,1));
         hs_max=h_max_vec(hmax_index,:);
     else
         hs_max=[rad(re_max),lil_h_max_re];
     end
     %and now we need to repeat this for finding the error in h
     val_min=1-(pi*rad.^2/a_p).*(var_dat(:,2)-var_dat(:,3));
     re_min=find(val_min > 0) ; im_min=find(val_min < 0);
     im_min_count=numel(im_min);
     val1_min=val_min(re_min);
     lil_h_min_re=rad(re_min)-sqrt(val1_min).*rad(re_min);
     if im_min_count > 0
         val2_min=val_min(im_min);
         lil_h_min_im=rad(im_min)+sqrt(abs(val2_min)).*rad(im_min);
         h_min_vec=[rad(re_min),lil_h_min_re;...
                    rad(im_min),lil_h_min_im];
         [hmin_srt,hmin_index]=sort(h_min_vec(:,1));
         hs_min=h_min_vec(hmin_index,:);
     else
         hs_min=[rad(re_min),lil_h_min_re];
     end
     %and now we calculate the error
     my_h=(hs_min(:,2)+hs_max(:,2))/2;
     del_h=hs_max(:,2)-my_h;
     error=del_h;
 end

end