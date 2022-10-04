clear all
close all
% format long g

file_path='F:\Chieco\Hyperuniformity\HU Jammed Packings from Carl Goodrich';
% tol_vec={'\tol1.0E-12'};
jam_type='\QuenchFromTinf';


%We have and exact expression for the average area or a particle when we
%are in the continuum limit. We compute the area here, knowing the ratio
%between the two species of particles is 1.4
x_min=0; x_max=1; y_min=0; y_max=1;
L_sys=x_max-x_min; sig_con=1.4;
A_sys=L_sys^2;
sys_specs=[x_min,y_min;x_max,y_max];
p_f_vec=[8400,8420];
for p_f=1:1:numel(p_f_vec)
    p_frac=p_f_vec(p_f);
    txt_phi=num2str(p_frac);
    fname=strcat(file_path,jam_type,'\particle positions',...
        '\0p',txt_phi,'xyrMats_200.txt');
    mat=dlmread(fname);
    r_start=1 ;r_fin=20;%numel(mat(1,:))/3;
    %first we read in each file
    for read=r_start:r_fin        
        if read==r_start
            N_tot=numel(mat(:,1));
            mat_stack=zeros(N_tot,3,r_fin);
            N_big=N_tot/2; N_small=N_tot/2;
        end
        mat_stack(:,:,read)=mat(:,3*(read-1)+1:3*(read-1)+3);
    end
    phi=p_frac/10^(ceil(log10(p_frac+1)));
    norm=phi
    %Now we compute diameter of a small particle, sig_s
    sig_s=sqrt(8*phi*A_sys/(N_tot*pi*(1+sig_con^2)));
    a_small=pi*(sig_s/2)^2;
    a_big=pi*(sig_con*sig_s/2)^2;
    %here is the the phi weighted average area of a particle, <a>
    phi_big=N_big*a_big/A_sys; phi_small=N_small*a_small/A_sys;
    avg_a=(phi_big*a_big+phi_small*a_small)/phi;
    a_vec=[zeros(N_tot/2,1)+a_small;zeros(N_tot/2,1)+a_big];
    avg_a_cubed=sum(a_vec.^4)/phi;
%     png_save=['D:\Chieco\My Papers\hyperuniformity soft discs\images\unedited images' jam_type ' phi_' txt_phi '.png'];
%     my_try=display_configs([mat_stack(:,1:2,1),[zeros(N_tot/2,1)+sig_s/2;zeros(N_tot/2,1)+(1.4*sig_s)/2]],...
%                                [0,0;1,1],png_save);
    c_rad_tenthousand=sqrt(L_sys^2/1E4);
    %we find how many windows and hat values here
    %         c_rads2=L_sys./(L_sys/c_rad_tenthousand:-1:2);
    %         c_rads1=logspace(log10(1E-2*sqrt(avg_a)),log10(c_rads2(1)-c_rads2(1)/2),15);
    %         c_rads3=logspace(log10(c_rads2(end)+c_rads2(end)/2),log10(0.99*L_sys),15);
    %         c_rads=[c_rads1,c_rads2,c_rads3];
    c_rads1=logspace(log10(0.5E-3*sqrt(avg_a)),log10(sqrt(avg_a)/2),15);
    c_rads2=logspace(log10(sqrt(avg_a)/2),log10(0.4*L_sys),30);
    c_rads=[c_rads1,c_rads2(2:end)];        
    N_windows_vec=zeros(1,numel(c_rads))+1E3;%round(L_sys^2./(c_rads.^2));
    %Now we test for hyperuniformity at a given box size L and keep all of
    %the raw data in a matrix
    i1_start=1; i1_fin=numel(c_rads);
    vec_write=zeros(i1_fin,5);
    mat_write_stack=zeros(i1_fin,5,r_fin);
    for i1=i1_start:i1_fin
        n_in=N_windows_vec(i1);
        L_in=c_rads(i1);
        n_phi_vec=r_fin*n_in;
        phi_vec=zeros(n_phi_vec,1);
        phi_sq_vec=phi_vec;
        for i2=r_start:r_fin
            phi_list=circle_window_periodic(sys_specs,[mat_stack(:,1:2,i2),a_vec],n_in,L_in,'Random');
            mat_start=(i2-1)*n_in+1; mat_fin=i2*n_in;
            phi_vec(mat_start:mat_fin,1)=phi_list(:,1);
            phi_sq_vec(mat_start:mat_fin,1)=phi_list.^2;
            var_matlab_loop=var(phi_list);
            mat_write_stack(i1,:,i2)=[L_in/sqrt(avg_a),avg_a,avg_a_cubed,n_in,var_matlab_loop/norm];
        end
        %Now we compute the 'corrected two pass variance' from Numerical
        %recipes. It is also called 'compensated variant' on wikipedia.
        %This is what Matlab uses to calculate the variance
        var_matlab=var(phi_vec);
        vec_write(i1,:)=[L_in/sqrt(avg_a),avg_a,avg_a_cubed,n_in,var_matlab/norm];
        text_path1=strcat(file_path,jam_type,'\real space averages\variance');
        filename=['\variance calc_from_N ' txt_phi '_avg.txt'];
        dlmwrite(strcat(text_path1,filename),vec_write(1:i1,:),'Delimiter',',','Precision','%1.15e','newline','pc');
    end
    text_path2=strcat(file_path,jam_type,'\loop data\variance');
    for i3=r_start:r_fin
        var_mat_write=mat_write_stack(:,:,i3);
        filename=['\variance calc_from_N' txt_phi ' loop_',num2str(i3),'.txt'];
        dlmwrite(strcat(text_path2,filename),var_mat_write,'Delimiter',',','Precision','%1.15e','newline','pc');
    end
    bark=0;
end
