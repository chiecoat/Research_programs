% clear all
% close all


person='HU Jammed Packings from Carl Goodrich';

%Weight by ys%
if strcmp(person,'HU Jammed Packings from Ning Xu') == 1
    titles={'10','20','30','40','100'};
    type='\quenchFromTinf\tol1.0E-12';
    phi_c_vec=[0.842,0.842,0.842,0.842,0.842];
    factor=20;
    pow_fix=1.8;
%     power_fit=@(a,xdata)(a(2)*((a(1)-xdata)+abs(a(1)-xdata)));
    fin_fit=@(a_out,phiC_m_phi)(a_out(1)*phiC_m_phi.^(1.8));
    phi_left=[0.828,0.828,0.828,0.828,0.828];
    phiC_min_phi_minimum=[0.0022,0.0022,0.0022,0.0022,0.0022];%zeros(1,5)+1E-3;%
    L_c_min=(0.00039223)/(323/43);
end

if strcmp(person,'HU Jammed Packings from Carl Goodrich') == 1
    type='quenchFromTinf';
    titles={'10','20','30','40'};
    if strcmp(type,'quenchFromTinf') == 1
        phi_c_vec=[0.8409,0.8409,0.8409,0.8409];
        factor=1;
        pow_fix=1.8;
%         power_fit=@(a,xdata)(a(2)*((a(1)-xdata)+abs(a(1)-xdata)));
        fin_fit=@(a_out,phiC_m_phi)(a_out(1)*phiC_m_phi.^(1.8));  
        phi_left=[0.804,0.804,0.804,0.804]; 
        phiC_min_phi_minimum=[1E-2,1E-2,1E-2,1E-2];
        L_c_min=0.00039223;
    else
        phi_c_vec=[0.8465,0.8465,0.8465,0.8465];
        factor=0.25;
        pow_fix=3/2;
%         power_fit=@(a,xdata)(a(2)*(a(1)-xdata)).*((a(1)-xdata)+abs(a(1)-xdata));
        fin_fit=@(a_out,phiC_m_phi)(a_out(1)*phiC_m_phi.^pow_fix);
        phi_left=[0.70,0.70,0.70,0.70];   
        phiC_min_phi_minimum=[0.02,0.015,0.023,0.02];
        L_c_min=0.00039223;
    end
end

pow_coeff=zeros(1,4);

for i1=1:numel(titles)
    file_path=strcat('F:\Chieco\Hyperuniformity\', ...
        person,'\', type,'\real space averages\normalized h difference');
    data_path=strcat('\h_min_he_norm L_root_a=',titles{i1},'.txt');
    filename=strcat(file_path,data_path);
    xy=csvread(filename,1,0);
    phi_c=phi_c_vec(i1);
    mat_fit=xy(xy(:,1)>=phi_left(i1),:);
    mat_in_x=mat_fit(mat_fit(:,2)>phiC_min_phi_minimum(i1),:);
    mat_in=mat_in_x(mat_in_x(:,3)>L_c_min,:);
    xs=mat_in(:,2);
    ys=mat_in(:,3);
    y_err=mat_in(:,4);
    y_weights=1./(y_err.^2);
    my_as=[factor];
    %compute the fit
    mdl=NonLinearModel.fit(xs,ys,fin_fit,my_as,'Weights',y_weights);
    coeffs_new=mdl.Coefficients.Estimate
    coeffs_error=mdl.Coefficients.SE;
    xstart=phi_left(i1); xfin=phi_c;
    xfit=logspace(-4,log10(xfin-xstart));
    y_plot=coeffs_new(1)*(xfit).^pow_fix;
%     if i1==1
%         loglog([1E-3,max(xs)],[4,4]*1E-4)
%         hold on
%     end
    errorbar(xs,ys,y_err,'o');
    set(gca,'XScale','log','YScale','log')
    hold on
    plot(xfit,y_plot)
    file_write=strcat('\power fits phiC=842 L_eq_',titles{i1},' root_a.txt');
    csvwrite(strcat(file_path,file_write),[xfit',y_plot'])
    bark=0;
end