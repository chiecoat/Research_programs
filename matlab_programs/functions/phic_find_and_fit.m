clear all
close all


person='HU Jammed Packings from Ning Xu';

%Weight by ys%
if strcmp(person,'HU Jammed Packings from Ning Xu') == 1
    titles={'30'};
    type='quenchFromTinf\tol1.0E-12';
    phi_c_vec=[0.842,0.842,0.842];
    factor=4;
    power_fit=@(a,xdata)(a(2)*((a(1)-xdata)+abs(a(1)-xdata)));
    fin_fit=@(a_out,phiC_m_phi)(a_out(1)*((phiC_m_phi)+abs(phiC_m_phi)));
    phi_left=0.78;
end

if strcmp(person,'HU Jammed Packings from Carl Goodrich') == 1
    type='quenchFromTinf';
    titles={'10','20','30'};
    if strcmp(type,'quenchFromTinf') == 1
        phi_c_vec=[0.842,0.842,0.842];
        factor=1;
        power_fit=@(a,xdata)(a(2)*((a(1)-xdata)+abs(a(1)-xdata)));
        fin_fit=@(a_out,phiC_m_phi)(a_out(1)*((phiC_m_phi)+abs(phiC_m_phi)));  
        phi_left=0.78; 
    else
        phi_c_vec=[0.85,0.85,0.85];
        factor=0.25;
        power_fit=@(a,xdata)(a(2)*(a(1)-xdata)).*((a(1)-xdata)+abs(a(1)-xdata));
        fin_fit=@(a_out,phiC_m_phi)(a_out(1)*phiC_m_phi).*((phiC_m_phi)+abs(phiC_m_phi));
        phi_left=0.40;    
    end
end

for i1=1:numel(titles)
    file_path=strcat('F:\Chieco\Hyperuniformity\', ...
        person,'\', type,'\real space averages\normalized h difference');
    data_path=strcat('\h_min_he_norm L_root_a=',titles{i1},'.txt');
    filename=strcat(file_path,data_path);
    xy=csvread(filename,1,0);
    bark=0;
    for i2 =1:numel(xy(:,1))-2
        if i2==1
            mat_to_plot=zeros(numel(xy(:,1))-2,3);
        end
        %Now we fit the data
         xs=xy(and(xy(:,1)>=phi_left,xy(:,1)<=1),1);
         ys=xy(and(xy(:,1)>=phi_left,xy(:,1)<=1),3);
         y_err=xy(and(xy(:,1)>=phi_left,xy(:,1)<=1),4);
%         xs=xy(i2:numel(xy(:,1)),1);
%         ys=xy(i2:numel(xy(:,1)),2);
%         y_err=xy(i2:numel(xy(:,1)),3);
        %we weight by the error
        y_weights=1./(y_err.^2);
        as=[phi_c_vec(i1),factor];
        %compute the fit
        mdl=NonLinearModel.fit(xs,ys,power_fit,as,'Weights',y_weights);
        coeffs=mdl.Coefficients.Estimate;
        err=mdl.Coefficients.SE;
        %plot the fit and compare to data
        xstart=min(xs); xfin=0.9;
        xfit=(xstart:(xfin-xstart)/500:xfin);
        if strcmp(type,'quenchFromTinf') == 1
            y_plot=coeffs(2)*((coeffs(1)-xfit)+abs(coeffs(1)-xfit));
        else
            y_plot=coeffs(2)*(coeffs(1)-xfit).*((coeffs(1)-xfit)+abs(coeffs(1)-xfit));
        end
        if i1==1
            errorbar(xy(:,1),xy(:,2),xy(:,3),'or')
            hold on
            plot(xfit,y_plot,'b')
            mat_to_plot(i2,:)=[min(xs),coeffs(1),err(1)];
            bark=0;
        end
        if i1==2
            errorbar(xy(:,1),xy(:,2),xy(:,3),'+r')
            hold on
            plot(xfit,y_plot,'g')
            mat_to_plot(i2,:)=[min(xs),coeffs(1),err(1)];
        end
    end
    bark=0;
    mat_to_fit=mat_to_plot(and(mat_to_plot(:,3)>0,mat_to_plot(:,1)>phi_left),:);
    phi_c_try=fitlm(mat_to_fit(:,1),mat_to_fit(:,2),'constant');
    phi_c=phi_c_try.Coefficients.Estimate;
    phi_c_err=fitlm(mat_to_fit(:,1),mat_to_fit(:,2),'linear','Weights',1./(mat_to_fit(:,3).^2));
    lin_stuff=phi_c_err.Coefficients.Estimate;
    hold off
    figure(i1)
    plot(mat_to_fit(:,1),mat_to_fit(:,2),'+b')
    hold on
    plot_me=[mat_to_fit(:,1),lin_stuff(1)+lin_stuff(2).*mat_to_fit(:,1)];  
    plot(plot_me(:,1),plot_me(:,2),'or')
    err=max(abs([max(plot_me(:,2))-phi_c,min(plot_me(:,2))-phi_c]));
    plot(plot_me(:,1),0.*plot_me(:,1)+phi_c,plot_me(:,1),0.*plot_me(:,1)+phi_c-err,plot_me(:,1),0.*plot_me(:,1)+phi_c+err,'g')
    bark=0;
    if strcmp(type,'quenchFromTinf') == 1
        mat_in=xy(and(xy(:,1)<=phi_c,xy(:,1)>=phi_left),:);
        xs=phi_c-mat_in(:,1);
        ys=mat_in(:,2);
        y_err=mat_in(:,3);
        %we weight by the error
        y_weights=1./(y_err.^2);%y_err;%
        my_as=[factor];
        %compute the fit
        mdl=NonLinearModel.fit(xs,ys,fin_fit,my_as,'Weights',y_weights);
        coeffs_new=mdl.Coefficients.Estimate;   
        xstart=min(xs); xfin=phi_c;
        xfit=logspace(-4,log10(max(xs(:,1))));
        y_plot=coeffs_new(1)*((xfit)+abs(xfit));
    else
        mat_in=xy(and(xy(:,1)<=phi_c,xy(:,1)>=phi_left),:);
        xs=phi_c-mat_in(:,1);
        ys=mat_in(:,2);
        y_err=mat_in(:,3);
        %we weight by the error
        y_weights=1./(y_err.^2);%y_err;%
        my_as=[factor];
        %compute the fit
        mdl=NonLinearModel.fit(xs,ys,fin_fit,my_as,'Weights',y_weights);
        coeffs_new=mdl.Coefficients.Estimate;       
        xstart=min(xs); xfin=phi_c;
        xfit=logspace(-4,log10(max(xs)));
        y_plot=coeffs_new(1)*(xfit).*((xfit)+abs((xfit)));
    end
%     file_write=strcat('\new power fits try L_eq_',titles{i1},' root_a.txt');
%     csvwrite(strcat(file_path,file_write),[xfit',y_plot'])  
%     [phi_c,err]
%     mdl
%     bark=0;
end

