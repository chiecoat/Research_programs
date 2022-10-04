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
function [lil_h_vec_write,error]=calculate_h_square(var_dat,a_p)

 big_l=var_dat(:,1);
%now we calculate h. Because of the nature of small window statistic there may be some L
%values where h might b imaginary. val will help compensate for that
 val=1-(big_l.^2/a_p).*var_dat(:,2);
 re=find(val > 0) ; im=find(val < 0);
 im_count=numel(im);
 val1=val(re);
 lil_h_re=big_l(re)/2-sqrt(val1).*big_l(re)/2;
 if im_count > 0 
   val2=val(im);
   lil_h_im=big_l(im)/2+sqrt(abs(val2)).*big_l(im)/2;
   lil_h_vec=[big_l(re),lil_h_re;...
              big_l(im),lil_h_im];
   [ell_sort,ell_index]=sort(lil_h_vec(:,1));
   lil_h_vec_write=lil_h_vec(ell_index,:);
 else
   lil_h_vec_write=[big_l(re),lil_h_re];
 end
 bark=0;

 %here we calculate the error in our hu length, if requested
 if nargout == 2
     %and now we need to repeat this for finding the error in h
     val_max=1-(big_l.^2/a_p).*(var_dat(:,2)+var_dat(:,3));
     re_max=find(val_max > 0) ; im_max=find(val_max < 0);
     im_max_count=numel(im_max);
     val1_max=val_max(re_max);
     lil_h_max_re=big_l(re_max)/2-sqrt(val1_max).*big_l(re_max)/2;
     if im_max_count > 0 
         val2_max=val_max(im_max);
         lil_h_max_im=big_l(im_max)/2+sqrt(abs(val2_max)).*big_l(im_max)/2;
         h_max_vec=[big_l(re_max),lil_h_max_re;...
                    big_l(im_max),lil_h_max_im];
         [hmax_sort,hmax_index]=sort(h_max_vec(:,1));
         hs_max=h_max_vec(hmax_index,:);
     else
         hs_max=[big_l(re_max),lil_h_max_re];
     end
     %and now we need to repeat this for finding the error in h
     val_min=1-(big_l.^2/a_p).*(var_dat(:,2)-var_dat(:,3));
     re_min=find(val_min > 0) ; im_min=find(val_min < 0);
     im_min_count=numel(im_min);
     val1_min=val_min(re_min);
     lil_h_min_re=big_l(re_min)/2-sqrt(val1_min).*big_l(re_min)/2;
     if im_min_count > 0
         val2_min=val_min(im_min);
         lil_h_min_im=big_l(im_min)/2+sqrt(abs(val2_min)).*big_l(im_min)/2;
         h_min_vec=[big_l(re_min),lil_h_min_re;...
                    big_l(im_min),lil_h_min_im];
         [hmin_srt,hmin_index]=sort(h_min_vec(:,1));
         hs_min=h_min_vec(hmin_index,:);
     else
         hs_min=[big_l(re_min),lil_h_min_re];
     end
     %and now we calculate the error
     my_h=(hs_min(:,2)+hs_max(:,2))/2;
     del_h=hs_max(:,2)-my_h;
     error=del_h;
 end

end
