% I think we’ll want to use a hybrid method where exactly-equal q values 
% for chi are averaged together at smallest q (crucial not to introduce 
% any systematic errors in this regime) and values falling in {q±dq} rings 
% are averaged at larger q (this won’t introduce noticeable systematic errors 
% and will beat down statistical uncertainty).
% 
% Could be done as follow:
% 
% — pick a dq<1 specifying rings [q-dq/2, q+dq/2) or [q, q+dq)
% 
% — Compute <q>, stdev[q], <chi>, stdev[chi] for points in each ring
% 
% — Record all four values only for nonempty rings, and use stdev values 
% as x- and y-error bars (at least for us to look at, but maybe for publication plots too?)
% 
% — What dq gives exactly-equal-q values method for first 5-10 points?
% 
% [1] rings could be logarithmically spaced:  [q_lo, q_lo(1+f) )
% 
% [2] If rate-limiting step is radially averaging, then do could do ensemble
% average by computing radial averaging on <FFT image> = (sum of FFT images)/Nensemble

function fft_radial_out=fft_radial_large_edit(im1,norm1,int_rad_max,dq,sqrt_avga)
% hold off
% Image size and normalization
rows=numel(im1(:,1));
cols=numel(im1(1,:));

% Fourier transform
F1 = fftshift(fft2(fftshift(full(im1))));
% F1 = fftshift(fft2(full(im1)));
clear im1
F = (abs(F1*norm1).^2);
% imshow(log(F),[ceil(min(min(log(F)))) ceil(max(max(log(F))))]); colormap(jet); colorbar
% clear F1

% Show spectrum (log)
% imagesc(log(F)); title('Fourier transform (abs log)'); 


freq=(-cols/2:cols/2-1);
[X, Y]=meshgrid(freq,freq); % get coordinates of the power spectrum image
[theta, rho]=cart2pol(X,Y); % equivalent in polar coordinates

n_max=rows/2-1; %max value we can use is rows/2-1

%%This is for every bin possible
% r_vec_all=unique(rho);
% 
% r_vec_all=r_vec_all(r_vec_all<n_max);
% chi_q_mat_all=zeros(numel(r_vec_all),4);
% 
% for r=1:numel(r_vec_all)
%     vals=F(find(rho==r_vec_all(r)));
%     chi_q_mat_all(r,1)=r_vec_all(r);
%     chi_q_mat_all(r,2)=std(r_vec_all(r));
%     chi_q_mat_all(r,3)=mean(vals); % average power values of all the same freq
%     chi_q_mat_all(r,4)=std(vals);
% %     if r_vec_all(r)>0
% %         loglog(zeros(numel(vals),1)+2*pi*sqrt_avga*chi_q_mat_all(r,1)/cols,vals,'.r','MarkerSize',12);
% %         hold on
% %     end
% end
% chi_all_plot=[2*pi*sqrt_avga*chi_q_mat_all(and(chi_q_mat_all(:,3)>0,chi_q_mat_all(:,3)<8),1:2)/cols,...
%               chi_q_mat_all(and(chi_q_mat_all(:,3)>0,chi_q_mat_all(:,3)<8),3:4)];
% % errorbar(chi_all_plot(:,1),chi_all_plot(:,3),chi_all_plot(:,4),chi_all_plot(:,4),chi_all_plot(:,2),chi_all_plot(:,2),...
% %          '-o','LineWidth',1.5,'Color',[153, 0, 0]/255,'CapSize',0.1)
% % set(gca,'YScale','log','XScale','log');
% loglog(chi_all_plot(:,1),chi_all_plot(:,3),'-o','LineWidth',1.5,'Color',[153, 0, 0]/255)
% hold on

% %Bins values between [r,r+dq) for constant dq
% chi_q_mat=zeros(n_max+1,2);
% %we know there are no values binned between 0 and 1, and we will have an
% %easier time getting the small q exactly right starting at 1.
% r_vec_step=[0,(1:dq:n_max+dq)];
% for r=1:numel(r_vec_step)-1
%     vals=F(find(and(rho>=r_vec_step(r),rho<r_vec_step(r+1))));
%     chi_q_mat(r,1)=r_vec_step(r);
%     chi_q_mat(r,2)=mean(vals); % average power values of all the same freq
% %     if and(r_vec(r)>0,numel(vals)>0)==1
% %         loglog(zeros(numel(vals),1)+2*pi*sqrt_avga*r_vec(r)/cols,vals,'+k','MarkerSize',4);
% %     end
% end
% loglog(2*pi*sqrt_avga*chi_q_mat(and(chi_q_mat(:,2)>0,chi_q_mat(:,2)<8),1)/cols,...
%      chi_q_mat(and(chi_q_mat(:,2)>0,chi_q_mat(:,2)<8),2),'-o','LineWidth',1.5,'Color',[155,155,155]/255);
 
%Bins values between [r-dq/2,r+dq/2) for constant dq

%we know there are no values binned between 0 and 1, and we will have an
%easier time getting the small q exactly right starting at 1.
% r_vec_step=[0,(1:dq:n_max+dq)];
% chi_q_mat_const=zeros(numel(r_vec_step),4);
% for r=1:numel(r_vec_step)-1
%     val_ind=find(and(rho>=r_vec_step(r)-dq/2,rho<r_vec_step(r)+dq/2));
%     chi_q_mat_const(r,1)=mean(rho(val_ind));
%     chi_q_mat_const(r,2)=std(rho(val_ind));
%     chi_q_mat_const(r,3)=mean(F(val_ind)); % average power values of all the same freq
%     chi_q_mat_const(r,4)=std(F(val_ind));
% %     if and(chi_q_mat_const(r,1)>0,numel(vals)>0)==1
% %         loglog(2*pi*sqrt_avga*rho(val_ind)/cols,F(val_ind),'.k','MarkerSize',12);
% %     end
% end
% chi_const_plot=[2*pi*sqrt_avga*chi_q_mat_const(and(chi_q_mat_const(:,3)>0,chi_q_mat_const(:,3)<8),1:2)/cols,...
%               chi_q_mat_const(and(chi_q_mat_const(:,3)>0,chi_q_mat_const(:,3)<8),3:4)];
% % errorbar(chi_const_plot(:,1),chi_const_plot(:,3),chi_const_plot(:,4),chi_const_plot(:,4),chi_const_plot(:,2),chi_const_plot(:,2),...
% %          '-o','LineWidth',1.5,'Color',[153,153,153]/255,'CapSize',0.1)
% loglog(chi_const_plot(:,1),chi_const_plot(:,3),'-o','LineWidth',1.5,'Color',[153,153,153]/255)

%Bins values between [r-dq/2,r+dq/2) for logarithmically spaced q values
%we know there are no values binned between 0 and 1, and we will have an
%easier time getting the small q exactly right starting at 1.
r=0;
count=0;
while r<n_max
    if count==0
        r=0;
        vals=F(find(and(rho>=r-dq/2,rho<r+dq/2)));
        chi_q_mat_log=[r,0,mean(vals),0];
        r_prime=1;
    elseif and(count>0,count<=int_rad_max)==1
        r=count;
        r_prime=r+1;
        val_log=find(and(rho>=r,rho<r_prime));
        %we want the same number for points for all runs regardless of if
        %we find values in the ring for any particular run. 
        if numel(val_log)>0
            r_max=max(rho(val_log)); r_min=min(rho(val_log));
            del_r=(r_max-r_min)/2;
            del_chi=std(F(val_log))/sqrt(numel(F(val_log)));
            chi_q_mat_log=[chi_q_mat_log;mean(rho(val_log)),del_r,mean(F(val_log)),del_chi];
        else
            chi_q_mat_log=[chi_q_mat_log;mean(rho(val_log)),std(rho(val_log)),mean(F(val_log)),std(F(val_log))];
        end
    else
        r=r_prime;
        r_prime=r*(1+dq);
        val_log=find(and(rho>=r,rho<r_prime));
        if numel(val_log)>0
            r_max=max(rho(val_log)); r_min=min(rho(val_log));
            del_r=(r_max-r_min)/2;
            del_chi=std(F(val_log))/sqrt(numel(F(val_log)));
            chi_q_mat_log=[chi_q_mat_log;mean(rho(val_log)),del_r,mean(F(val_log)),del_chi];
        else
            chi_q_mat_log=[chi_q_mat_log;mean(rho(val_log)),std(rho(val_log)),mean(F(val_log)),std(F(val_log))];
        end
    end
    count=count+1;
%     if and(chi_q_mat_log(count)>0,numel(vals)>0)==1
%         loglog(2*pi*sqrt_avga*rho(val_log)/cols,F(val_log),'.b','MarkerSize',12);
%     end
end
% chi_log_plot=[2*pi*sqrt_avga*chi_q_mat_log(and(chi_q_mat_log(:,3)>0,chi_q_mat_log(:,3)<8),1:2)/cols,...
%               chi_q_mat_log(and(chi_q_mat_log(:,3)>0,chi_q_mat_log(:,3)<8),3:4)];
% errorbar(chi_log_plot(:,1),chi_log_plot(:,3),chi_log_plot(:,4),chi_log_plot(:,4),chi_log_plot(:,2),chi_log_plot(:,2),...
%          '-o','LineWidth',1.5,'Color',[0,204,204]/255,'CapSize',0.1)
% loglog(chi_log_plot(:,1),chi_log_plot(:,3),'-o','LineWidth',1.5,'Color',[0,204,204]/255)
 

fft_radial_out=[2*pi*sqrt_avga*chi_q_mat_log(:,1:2)/cols,chi_q_mat_log(:,3:4)];