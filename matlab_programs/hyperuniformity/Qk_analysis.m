% +
% NAME: Qk_analysis
%
% PURPOSE:
%     This program checks that Qk is healthy. If Qk~=0 then we are good, it
%     should be incredibly close to zeros with only about roundoff error as
%     the difference from 0. If this is the case, and Dk is much larger
%     than zero, then our data is healthy. Then we can bin and histogram the 
%     which we will fit a histogram to.
%
% CATEGORY:
%     Hyperuniformity
%
% CALLING SEQUENCE:
%    Qk_analysis
%
% INPUTS: 
%    None but we'll need access to the files that the Qk code 
%
% OPTIONAL INPUTS: (none)
%
% KEYWORD PARAMETERS: (none)
%
% OUTPUTS: Will histogram Qk and fit a Gaussian to it. Also we will find
% the fraction of particles in Qk that participat below in
% underpacked/overpacked regions.
% 
% SIDE EFFECTS: (none)
%
% MODIFICATION HISTORY:
%    written by: A. Chieco, UPenn, September 2016
%
%-
filenames={'voro_phi850_1 divergence', 'step1e4_confire_85 divergence',...
           'step1e4_confire_87 divergence','step1e6_confire_87 divergence'};

for i1=1:numel(filenames)
    if i1==1
        Qk_path='D:\Chieco\Hyperuniformity\HU Jammed Packings from Ning Xu\quenchN1E6\Q_k_folder';
    else
        Qk_path='D:\Chieco\Hyperuniformity\HU Jammed Packings from Ning Xu\quenchVaryRate\Q_k analysis';
    end
    Qk_mat=load([Qk_path '\' filenames{i1} '.mat']);
    
    Qk=Qk_mat.QkStatic;
    Qk=Qk(abs(Qk)<=1);
    
    hist=histogram(Qk,200,'Normalization','pdf');
    set(gca,'YScale','log')
    ylim([1E-5,100]); xlim([-0.6,0.3])
    %First we have to find which particles are large and small
    %This is done by histogramming the "mass" column of the data
    edges=hist.BinEdges;
    xs=edges(1:numel(edges)-1);
    val=hist.Values;
    g_coeffs=[0.1];
    
    g_fit_xs=xs(xs<=0.3);
    g_fit_ys=val(xs<=0.3);
    my_gauss=@(coeffs,qs)(1/(sqrt(2*pi)*coeffs(1)))*exp(-((qs/(sqrt(2)*coeffs(1))).^2));
    my_gauss_coeffs=nlinfit(g_fit_xs,g_fit_ys,my_gauss,g_coeffs)
    plot_points_gauss=(1/(sqrt(2*pi)*my_gauss_coeffs(1)))*exp(-((xs/(sqrt(2)*my_gauss_coeffs(1))).^2));
    hold on
    plot(xs,plot_points_gauss,'r','LineWidth',2)
    
%     g_plus_xs=xs(xs>0.0);
%     g_plus_ys=val(xs>0.0);
%     my_plus_gauss=fit(g_plus_xs',g_plus_ys','gauss1');
%     mp_g_a1=my_plus_gauss.a1;
%     mp_g_b1=my_plus_gauss.b1;
%     mp_g_c1=my_plus_gauss.c1;
%     plot_plus_gauss=mp_g_a1*exp(-(((xs-mp_g_b1)/mp_g_c1).^2));
%     plot(xs,plot_plus_gauss,'g','LineWidth',2)
    
    %We fit for the parameters of the Guassian so we they are hard coded into
    %our "full" Qk analysis    
    x_pos=xs(val>0);
    val_pos=val(val>0);
    x_in=x_pos(x_pos<-0.1);
    val_in=val_pos(x_pos<-0.1);
    w_coeffs=[4E3,.024];
    x_pow=2;
    
    weibull=@(coeffs,qs)coeffs(1)*(((qs-abs(qs))/2).^x_pow).*exp(-(abs(qs)/coeffs(2)));
    options=optimset('MaxIter',1E3);
    my_weib_coeffs=nlinfit(x_in,val_in,weibull,w_coeffs,options,'Weights',1./val_in.^2)
    plot_points_weibull=my_weib_coeffs(1)*(((xs-abs(xs))/2).^x_pow).*exp(-(abs(xs)/my_weib_coeffs(2)));
    plot(xs,plot_points_weibull,'g','LineWidth',2)
%     Qk_sort=sort(-Qk(Qk<0));
%     my_try=wblfit(Qk_sort)
%     weib=(my_try(2)/(my_try(1)^my_try(2)))*(Qk_sort.^(my_try(2)-1)).*exp(-(Qk_sort/my_try(1)).^my_try(2));
% %     plot(-Qk_sort,weib,'LineWidth',2)
%     
%     my_weibull=(my_try(2)/(my_try(1)^my_try(2)))*((-xs(xs<0)).^(my_try(2)-1)).*exp(-((-xs(xs<0))/my_try(1)).^my_try(2));
%     plot_points_weibull=zeros(1,numel(plot_points_gauss),1);
%     plot_points_weibull(1:numel(my_weibull))=my_weibull;
%     plot(xs,plot_points_weibull,'LineWidth',2)
    
    
    Qk_fit_full=@(coeffs,qs)coeffs(1)*plot_points_gauss(val>0)+(1-coeffs(1))*plot_points_weibull(val>0);
    Qk_coeffs_in=[1-0.0052];
    x_pos=xs(val>0);
    val_pos=val(val>0);
    Qk_coeffs=nlinfit(x_pos,val_pos,Qk_fit_full,Qk_coeffs_in,'Weights',1./val_pos.^2)
    plot_points_Qk=Qk_coeffs(1)*plot_points_gauss+(1-Qk_coeffs(1))*plot_points_weibull;
    plot(xs,plot_points_Qk,'k','LineWidth',2)
    
    bark=0;
    close Figure 1
end

