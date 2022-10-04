clear all
close all
% format long g

file_path='D:\Chieco\Hyperuniformity\HU Jammed Packings from Ning Xu';
jam_type='\quenchVaryRate';
tol_vec={'\tol1.0E-12\N=1E6'};
num_in=1E6;
num_in_txt='1E6';


%We have and exact expression for the average area or a particle when we
%are in the continuum limit. We compute the area here, knowing the ratio
%between the two species of particles is 1.4
x_min=-0.5; x_max=0.5; y_min=-0.5; y_max=0.5;
L_sys=x_max-x_min; sig_con=1.4;
A_sys=L_sys^2;
sys_specs=[x_min,y_min;x_max,y_max];
quench_rate={'step_0'};
p_f_vec=[86];

%This is the diameter we are aiming for for the small particles.
ds_scale=15;
a_small=pi*(ds_scale/2)^2;
imi_xs_small=meshgrid((-floor(ds_scale/2):floor(ds_scale/2)));
imi_ys_small=imi_xs_small';
[theta_small,rho_small]=cart2pol(imi_xs_small,imi_ys_small);
imi_small=zeros(numel(rho_small(:,1)),numel(rho_small(1,:)));
imi_small(rho_small<=ds_scale/2)=1;
%This is just sig_con*ds_scale but I don't want to leave any chance for
%rounding errors
dl_scale=21;
a_large=pi*(dl_scale/2)^2;
imi_xs_large=meshgrid((-floor(dl_scale/2):floor(dl_scale/2)));
imi_ys_large=imi_xs_large';
[theta_large,rho_large]=cart2pol(imi_xs_large,imi_ys_large);
imi_large=zeros(numel(rho_large(:,1)),numel(rho_large(1,:)));
imi_large(rho_large<=dl_scale/2)=1;

for p_f=1:numel(quench_rate)
    r_start=1; r_fin=1;
    %first we read in each file
    for t1=1:numel(p_f_vec)
        tol=tol_vec{1};
        p_frac=p_f_vec(t1);
        txt_phi=num2str(p_f_vec(t1));
        for read=r_start:r_fin
            fname=[file_path jam_type tol  '\particle positions\' quench_rate{p_f}...
                '\confire_0.' txt_phi '_' num2str(read)];
            mat=dlmread(fname,'\t',[1 0 num_in 1]);
            if read==r_start
                N_tot=numel(mat(:,1));
                mat_stack=zeros(N_tot,2,r_fin);
                N_big=N_tot/2; N_small=N_tot/2;
            end
            mat_stack(:,:,read)=mat;
        end
        phi=p_frac/10^(ceil(log10(p_frac+1)));
        norm=phi;
        %Now we compute diameter of a small particle, sig_s
        A_sys=1.0;
        sig_s=sqrt(8*phi*A_sys/(N_tot*pi*(1+sig_con^2)));
        sz_scale=ds_scale/sig_s;
        L_sys=ceil(sz_scale); A_sys=sz_scale^2;
        L_sys
        %We need two vectors for each particle population.
        mat_in=[sz_scale*[double(mat(:,1)-x_min),double(mat(:,2)-y_min)],...
            [zeros(N_small,1)+ds_scale/2;zeros(N_big,1)+dl_scale/2]];
        mat_in=[ceil(mat_in(:,1)),ceil(mat_in(:,2)),floor(mat_in(:,3))];
        imi=zeros(L_sys,L_sys);
        imi_name=[file_path jam_type tol '\phi85_fill.png'];
% % % %         for place_small=1:N_small
% % % %             xyr_vec=mat_in(place_small,:);
% % % %             %Now we have to check if our particle lies at all over the x
% % % %             %boundary
% % % %             x_left=xyr_vec(1)-xyr_vec(3); x_right=xyr_vec(1)+xyr_vec(3);
% % % %             y_bot=xyr_vec(2)-xyr_vec(3); y_top=xyr_vec(2)+xyr_vec(3);
% % % %             if and(x_left>0,x_right<=L_sys)==1
% % % %                 %Then we have to check the three y boundary conditions, the
% % % %                 %particle is either in the bulk or it has area lying above
% % % %                 %or below the boundary
% % % %                 if and(y_bot>=1,y_top<=L_sys)==1
% % % %                     imi(y_bot:y_top,x_left:x_right)=imi(y_bot:y_top,x_left:x_right)+imi_small;
% % % %                 end
% % % %                 if and(y_bot>=1,y_top>L_sys)==1
% % % %                     imi(y_bot:L_sys,x_left:x_right)=imi(y_bot:L_sys,x_left:x_right)+imi_small(1:L_sys-y_bot+1,:);
% % % %                     imi(1:y_top-L_sys,x_left:x_right)=imi(1:y_top-L_sys,x_left:x_right)+imi_small(end-(y_top-L_sys)+1:end,:);
% % % %                 end
% % % %                 if and(y_bot<1,y_top<=L_sys)==1
% % % %                     imi(1:y_top,x_left:x_right)=imi(1:y_top,x_left:x_right)+imi_small(end-y_top+1:end,:);
% % % %                     imi(L_sys-abs(y_bot):L_sys,x_left:x_right)=imi(L_sys-abs(y_bot):L_sys,x_left:x_right)+imi_small(1:abs(y_bot)+1,:);
% % % %                 end
% % % %             else
% % % %                 %Now we check for particles that lie over the x edge or the
% % % %                 %system corner. We first differentiate the left side from
% % % %                 %the right side of the system
% % % %                 %Here particles lie over the left edge
% % % %                 if and(x_left<1,x_right<=L_sys)==1
% % % %                     %Now we check if the particles lies over the left edge 
% % % %                     %but in the bulk
% % % %                     if and(y_bot>=1,y_top<=L_sys)==1
% % % %                         imi(y_bot:y_top,1:x_right)=imi(y_bot:y_top,1:x_right)+imi_small(:,end-x_right+1:end);
% % % %                         imi(y_bot:y_top,L_sys-abs(x_left):L_sys)=imi(y_bot:y_top,L_sys-abs(x_left):L_sys)+imi_small(:,1:abs(x_left)+1);
% % % %                     end
% % % %                     %Now we check to see if we are in the top left corner
% % % %                     if and(y_bot<1,y_top<=L_sys)==1
% % % %                         imi(1:y_top,1:x_right)=imi(1:y_top,1:x_right)+imi_small(end-y_top+1:end,end-x_right+1:end);
% % % %                         imi(1:y_top,L_sys-abs(x_left):L_sys)=imi(1:y_top,L_sys-abs(x_left):L_sys)+imi_small(end-y_top+1:end,1:abs(x_left)+1);
% % % %                         imi(L_sys-abs(y_bot):L_sys,1:x_right)=imi(L_sys-abs(y_bot):L_sys,1:x_right)+imi_small(1:abs(y_bot)+1,end-x_right+1:end);
% % % %                         imi(L_sys-abs(y_bot):L_sys,L_sys-abs(x_left):L_sys)=imi(L_sys-abs(y_bot):L_sys,L_sys-abs(x_left):L_sys)+imi_small(1:abs(y_bot)+1,1:abs(x_left)+1);
% % % %                     end
% % % %                     %Now we check to see if we are in the bottom left corner
% % % %                     if and(y_bot>=1,y_top>L_sys)==1
% % % %                         imi(y_bot:L_sys,1:x_right)=imi(y_bot:L_sys,1:x_right)+imi_small(1:L_sys-y_bot+1,end-x_right+1:end);
% % % %                         imi(y_bot:L_sys,L_sys-abs(x_left):L_sys)=imi(y_bot:L_sys,L_sys-abs(x_left):L_sys)+imi_small(1:L_sys-y_bot+1:end,1:abs(x_left)+1);
% % % %                         imi(1:y_top-L_sys,1:x_right)=imi(L_sys-abs(y_bot):L_sys,1:x_right)+imi_small(end-(y_top-L_sys)+1:end,end-x_right+1:end);
% % % %                         imi(1:y_top-L_sys,L_sys-abs(x_left):L_sys)=imi(L_sys-abs(y_bot):L_sys,L_sys-abs(x_left):L_sys)+imi_small(end-(y_top-L_sys)+1:end,1:abs(x_left)+1);
% % % %                     end
% % % %                 end
% % % %                 %Here particles lie over the right edge
% % % %                 if and(x_left>=1,x_right>L_sys)==1
% % % %                     %Now we check if the particles lies over the left edge 
% % % %                     %but in the bulk                  
% % % %                     if and(y_bot>=1,y_top<=L_sys)==1
% % % %                         imi(y_bot:y_top,x_left:L_sys)=imi(y_bot:y_top,x_left:L_sys)+imi_small(:,1:L_sys-x_left+1);
% % % %                         imi(y_bot:y_top,1:x_right-L_sys)=imi(y_bot:y_top,1:x_right-L_sys)+imi_small(:,end-(x_right-L_sys)+1:end);
% % % %                     end
% % % %                     %Now we check to see if we are in the top right corner
% % % %                     if and(y_bot<1,y_top<=L_sys)==1
% % % %                         imi(1:y_top,x_left:L_sys)=imi(1:y_top,x_left:L_sys)+imi_small(end-y_top+1:end,1:L_sys-x_left+1);
% % % %                         imi(1:y_top,1:x_right-L_sys)=imi(1:y_top,1:x_right-L_sys)+imi_small(end-y_top+1:end,end-(x_right-L_sys)+1:end);
% % % %                         imi(L_sys-abs(y_bot):L_sys,x_left:L_sys)=imi(L_sys-abs(y_bot):L_sys,x_left:L_sys)+imi_small(1:abs(y_bot)+1,1:L_sys-x_left+1);
% % % %                         imi(L_sys-abs(y_bot):L_sys,L_sys-abs(x_left):L_sys)=imi(L_sys-abs(y_bot):L_sys,L_sys-abs(x_left):L_sys)+imi_small(1:abs(y_bot)+1,end-(x_right-L_sys)+1:end);
% % % %                     end
% % % % %                   %Now we check to see if we are in the bottom right corner
% % % %                     if and(y_bot>=1,y_top>L_sys)==1
% % % %                         imi(y_bot:L_sys,x_left:L_sys)=imi(y_bot:L_sys,x_left:L_sys)+imi_small(1:L_sys-y_bot+1,1:L_sys-x_left+1);
% % % %                         imi(y_bot:L_sys,1:x_right-L_sys)=imi(y_bot:L_sys,1:x_right-L_sys)+imi_small(1:L_sys-y_bot+1:end,end-(x_right-L_sys)+1:end);
% % % %                         imi(1:y_top-L_sys,x_left:L_sys)=imi(1:y_top-L_sys,x_left:L_sys)+imi_small(end-(y_top-L_sys)+1:end,1:L_sys-x_left+1);
% % % %                         imi(1:y_top-L_sys,1:x_right-L_sys)=imi(L_sys-abs(y_bot):L_sys,1:x_right-L_sys)+imi_small(end-(y_top-L_sys)+1:end,end-(x_right-L_sys)+1:end);
% % % %                     end
% % % %                 end
% % % %             end
% % % %         end
% % % %         for place_large=N_big+1:N_tot
% % % %             xyr_vec=mat_in(place_large,:);
% % % %             %Now we have to check if our particle lies at all over the x
% % % %             %boundary
% % % %             x_left=xyr_vec(1)-xyr_vec(3); x_right=xyr_vec(1)+xyr_vec(3);
% % % %             y_bot=xyr_vec(2)-xyr_vec(3); y_top=xyr_vec(2)+xyr_vec(3);
% % % %             if and(x_left>0,x_right<=L_sys)==1
% % % %                 %Then we have to check the three y boundary conditions, the
% % % %                 %particle is either in the bulk or it has area lying above
% % % %                 %or below the boundary
% % % %                 if and(y_bot>=1,y_top<=L_sys)==1
% % % %                     imi(y_bot:y_top,x_left:x_right)=imi(y_bot:y_top,x_left:x_right)+imi_large;                    
% % % %                 end
% % % %                 if and(y_bot>=1,y_top>L_sys)==1
% % % %                     imi(y_bot:L_sys,x_left:x_right)=imi(y_bot:L_sys,x_left:x_right)+imi_large(1:L_sys-y_bot+1,:);
% % % %                     imi(1:y_top-L_sys,x_left:x_right)=imi(1:y_top-L_sys,x_left:x_right)+imi_large(end-(y_top-L_sys)+1:end,:);
% % % %                 end
% % % %                 if and(y_bot<1,y_top<L_sys)==1
% % % %                     imi(1:y_top,x_left:x_right)=imi(1:y_top,x_left:x_right)+imi_large(end-y_top+1:end,:);
% % % %                     imi(L_sys-abs(y_bot):L_sys,x_left:x_right)=imi(L_sys-abs(y_bot):L_sys,x_left:x_right)+imi_large(1:abs(y_bot)+1,:);
% % % %                 end                
% % % %             else
% % % %                 %Now we check for particles that lie over the x edge or the
% % % %                 %system corner. We first differentiate the left side from
% % % %                 %the right side of the system
% % % %                 %Here particles lie over the left edge
% % % %                 if and(x_left<1,x_right<=L_sys)==1
% % % %                     %Now we check if the particles lies over the left edge 
% % % %                     %but in the bulk
% % % %                     if and(y_bot>=1,y_top<=L_sys)==1
% % % %                         imi(y_bot:y_top,1:x_right)=imi(y_bot:y_top,1:x_right)+imi_large(:,end-x_right+1:end);
% % % %                         imi(y_bot:y_top,L_sys-abs(x_left):L_sys)=imi(y_bot:y_top,L_sys-abs(x_left):L_sys)+imi_large(:,1:abs(x_left)+1);
% % % %                     end
% % % %                     %Now we check to see if we are in the top left corner
% % % %                     if and(y_bot<1,y_top<=L_sys)==1
% % % %                         imi(1:y_top,1:x_right)=imi(1:y_top,1:x_right)+imi_large(end-y_top+1:end,end-x_right+1:end);
% % % %                         imi(1:y_top,L_sys-abs(x_left):L_sys)=imi(1:y_top,L_sys-abs(x_left):L_sys)+imi_large(end-y_top+1:end,1:abs(x_left)+1);
% % % %                         imi(L_sys-abs(y_bot):L_sys,1:x_right)=imi(L_sys-abs(y_bot):L_sys,1:x_right)+imi_large(1:abs(y_bot)+1,end-x_right+1:end);
% % % %                         imi(L_sys-abs(y_bot):L_sys,L_sys-abs(x_left):L_sys)=imi(L_sys-abs(y_bot):L_sys,L_sys-abs(x_left):L_sys)+imi_large(1:abs(y_bot)+1,1:abs(x_left)+1);
% % % %                     end
% % % %                     %Now we check to see if we are in the bottom right corner
% % % %                     if and(y_bot>=1,y_top>L_sys)==1
% % % %                         imi(y_bot:L_sys,1:x_right)=imi(y_bot:L_sys,1:x_right)+imi_large(1:L_sys-y_bot+1,end-x_right+1:end);
% % % %                         imi(y_bot:L_sys,L_sys-abs(x_left):L_sys)=imi(y_bot:L_sys,L_sys-abs(x_left):L_sys)+imi_large(1:L_sys-y_bot+1:end,1:abs(x_left)+1);
% % % %                         imi(1:y_top-L_sys,1:x_right)=imi(L_sys-abs(y_bot):L_sys,1:x_right)+imi_large(end-(y_top-L_sys)+1,end-x_right+1:end);
% % % %                         imi(1:y_top-L_sys,L_sys-abs(x_left):L_sys)=imi(L_sys-abs(y_bot):L_sys,L_sys-abs(x_left):L_sys)+imi_large(end-(y_top-L_sys)+1,1:abs(x_left)+1);
% % % %                     end
% % % %                 end
% % % %                 %Here particles lie over the right edge
% % % %                 if and(x_left>=1,x_right>L_sys)==1
% % % %                     %Now we check if the particles lies over the right edge 
% % % %                     %but in the bulk                  
% % % %                     if and(y_bot>=1,y_top<=L_sys)==1
% % % %                         imi(y_bot:y_top,x_left:L_sys)=imi(y_bot:y_top,x_left:L_sys)+imi_large(:,1:L_sys-x_left+1);
% % % %                         imi(y_bot:y_top,1:x_right-L_sys)=imi(y_bot:y_top,1:x_right-L_sys)+imi_large(:,end-(x_right-L_sys)+1:end);
% % % %                     end
% % % %                     %Now we check to see if we are in the top right corner
% % % %                     if and(y_bot<1,y_top<=L_sys)==1
% % % %                         imi(1:y_top,x_left:L_sys)=imi(1:y_top,x_left:L_sys)+imi_large(end-y_top+1:end,1:L_sys-x_left+1);
% % % %                         imi(1:y_top,1:x_right-L_sys)=imi(1:y_top,1:x_right-L_sys)+imi_large(end-y_top+1:end,end-(x_right-L_sys)+1:end);
% % % %                         imi(L_sys-abs(y_bot):L_sys,x_left:L_sys)=imi(L_sys-abs(y_bot):L_sys,x_left:L_sys)+imi_large(1:abs(y_bot)+1,1:L_sys-x_left+1);
% % % %                         imi(L_sys-abs(y_bot):L_sys,L_sys-abs(x_left):L_sys)=imi(L_sys-abs(y_bot):L_sys,L_sys-abs(x_left):L_sys)+imi_large(1:abs(y_bot)+1,end-(x_right-L_sys)+1:end);
% % % %                     end
% % % %                   %Now we check to see if we are in the bottom right corner
% % % %                     if and(y_bot>=1,y_top>L_sys)==1
% % % %                         imi(y_bot:L_sys,x_left:L_sys)=imi(y_bot:L_sys,x_left:L_sys)+imi_large(1:L_sys-y_bot+1,1:L_sys-x_left+1);
% % % %                         imi(y_bot:L_sys,1:x_right-L_sys)=imi(y_bot:L_sys,1:x_right-L_sys)+imi_large(1:L_sys-y_bot+1:end,end-(x_right-L_sys)+1:end);
% % % %                         imi(1:y_top-L_sys,x_left:L_sys)=imi(1:y_top-L_sys,x_left:L_sys)+imi_large(end-(y_top-L_sys)+1:end,1:L_sys-x_left+1);
% % % %                         imi(1:y_top-L_sys,1:x_right-L_sys)=imi(L_sys-abs(y_bot):L_sys,1:x_right-L_sys)+imi_large(end-(y_top-L_sys)+1:end,end-(x_right-L_sys)+1:end);
% % % %                     end
% % % %                 end
% % % %             end
% % % %         end 
        a_small=pi*(ds_scale/2)^2;
        a_big=pi*(dl_scale/2)^2;
        %here is the the phi weighted average area of a particle, <a>
        phi_big=N_big*a_big/A_sys; phi_small=N_small*a_small/A_sys;
        avg_a=(phi_big*a_big+phi_small*a_small)/phi;
        clear imi
        a_vec=[zeros(N_tot/2,1)+a_small;zeros(N_tot/2,1)+a_big];
        avg_a_cubed=sum(a_vec.^4)/phi;
        c_rad_tenthousand=sqrt(L_sys^2/1E4);
        %we find how many windows and hat values here
        c_rads1=(1:round(sqrt(avg_a)));
        c_rads2=logspace(log10(round(sqrt(avg_a))),log10(round(0.5*L_sys)),131-round(sqrt(avg_a)));
        c_rads=unique(round([c_rads1,c_rads2(2:end)]));        
        N_windows_vec=zeros(1,numel(c_rads))+1E4;        
        %Now we test for hyperuniformity at a given box size L and keep all of
        %the raw data in a matrix
        i1_start=1; i1_fin=numel(c_rads);
        vec_write=zeros(i1_fin,5);
        mat_write_stack=zeros(i1_fin,5,r_fin);
        clear imi
        for i1=1:i1_fin
            n_in=N_windows_vec(i1);
            L_in=c_rads(i1);
            n_phi_vec=r_fin*n_in;
            phi_vec=zeros(n_phi_vec,1);
            phi_sq_vec=phi_vec;
            for i2=r_start:r_fin
                phi_list=square_window_pixel_periodic(imi,n_in,L_in,'Random');
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
            text_path1=strcat(file_path,jam_type,tol,'\real space averages\variance');
            filename=['\' quench_rate{p_f} ' variance square_window_pixel imi_xl ' txt_phi '_avg.txt'];
            dlmwrite(strcat(text_path1,filename),vec_write(1:i1,:),'Delimiter',',','Precision','%1.15e','newline','pc');
        end
        text_path2=strcat(file_path,jam_type,tol,'\loop data\variance');
        for i3=r_start:r_fin
            var_mat_write=mat_write_stack(:,:,i3);
            filename=['\' quench_rate{p_f} ' variance square_window_pixel imi_xl ' txt_phi ' loop_',num2str(i3),'.txt'];
            dlmwrite(strcat(text_path2,filename),var_mat_write,'Delimiter',',','Precision','%1.15e','newline','pc');
        end
    end
end


