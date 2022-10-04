% clear all
% close all
% format long g

file_path='E:\Chieco\Hyperuniformity\HU foams';
data_type={'arc','vert'};       
imi_nums={'125','250','500'} ;
width_vec=1;
for dat=1:numel(data_type)
    for imi=1:numel(imi_nums)  
        num=imi_nums{imi};
        %we want to have one lengthscale that we normalize all the data by.
        %This will be th average bubble size
        foam_norm=dlmread(['E:\Chieco\Hyperuniformity\HU foams'...
                  '\patterns_to_analyze\bub_stats_mat imi_' num '.txt']);
        area_norm=sum(foam_norm(:,3).^2)/sum(foam_norm(:,3));
        bub_length_norm=sqrt(area_norm);
        %The system size is determined by the minimum and maximum positoins
        %of the vertices
        foam_sys=dlmread(['E:\Chieco\Hyperuniformity\HU foams\patterns_to_analyze'...
                           '\vert_stats_mat imi_' num  '.txt']);
        x_min=min(foam_sys(:,1)); x_max=max(foam_sys(:,1));
        y_min=min(foam_sys(:,2)); y_max=max(foam_sys(:,2));    
        %first we read in each file
        fname=[file_path '\patterns_to_analyze\' data_type{dat} '_stats_mat imi_' num '.txt'];        
        mat=dlmread(fname);    
        if strcmp(data_type{dat},'arc')==1
            len_pad=sum(mat(:,3).^2)/sum(mat(:,3));
            mat_keep_tot=mat(and(and(mat(:,1)>x_min+len_pad,mat(:,1)<x_max-len_pad),...
                and(mat(:,2)>y_min+len_pad,mat(:,2)<y_max-len_pad)),:);
            L_sys_x=x_max-x_min-2*len_pad;
            L_sys_y=y_max-y_min-2*len_pad;
            A_sys=L_sys_x*L_sys_y;
            L_sys=min([L_sys_x,L_sys_y]);
            sys_specs=[x_min,y_min;x_max,y_max];
        else
            mat_keep_tot=mat(and(and(mat(:,1)>x_min,mat(:,1)<x_max),...
                and(mat(:,2)>y_min,mat(:,2)<y_max)),:);
            L_sys_x=x_max-x_min;
            L_sys_y=y_max-y_min;
            A_sys=L_sys_x*L_sys_y;
            L_sys=min([L_sys_x,L_sys_y]);
            sys_specs=[x_min,y_min;x_max,y_max];
        end
        for w_count=1:numel(width_vec)
            width=width_vec(w_count);
            for runs=1:2%numel(mat(1,3:end))+1
                if runs==1
                    areas=mat_keep_tot(:,3);
                    type='points_weighted';
                    length_norm=sum(areas.^2)/sum(areas);%sqrt(avg_a);%
%                     a_vec=areas*width;
                    a_vec=areas*width;
                    phi=sum(a_vec)/A_sys;
                    norm=phi;
                elseif runs==numel(mat(1,3:end))+1
                    areas=mat_keep_tot(:,3)*0+1;
                    type='points_mono';
                    length_norm=sqrt(A_sys/numel(areas));
                    a_vec=areas;
                    phi=sum(a_vec)/A_sys;
                    norm=1;
                else
                    continue
                    %areas=mat(:,end);
                    %type='centers_coord';
                end                
                avg_a=sum(a_vec.^2)/sum(a_vec);
                avg_a_cubed=sum(a_vec.^4)/sum(a_vec);
                %we find how many windows and hat values here
                c_rads1=logspace(log10(0.5E-2*length_norm),log10(length_norm/2),25);
                c_rads2=logspace(log10(length_norm/2),log10(0.4*L_sys),5);
                c_rads=[c_rads1,c_rads2(2:end)];
                N_windows_vec=zeros(1,numel(c_rads))+1E5;%round(L_sys^2./(c_rads.^2));
                %Now we test for hyperuniformity at a given box size L and keep all of
                %the raw data in a matrix
                i1_start=1; i1_fin=numel(c_rads);
                vec_write=zeros(i1_fin,10);
                mat_write_stack=zeros(i1_fin,10);
                r_start=1; r_fin=1;
                for i1=1:i1_fin
                    n_in=N_windows_vec(i1);
                    L_in=c_rads(i1);
                    n_phi_vec=r_fin*n_in;
                    phi_vec=zeros(n_phi_vec,1);
                    phi_sq_vec=phi_vec;
                    for i2=r_start:r_fin
                        if runs==1
                            var_list=circle_window_centers(sys_specs,[mat_keep_tot(:,1:2),a_vec],n_in,L_in);
                        else
                            var_list=circle_window_points(sys_specs,mat_keep_tot(:,1:2),n_in,L_in); 
                        end
                        mat_start=(i2-1)*n_in+1; mat_fin=i2*n_in;
                        phi_vec(mat_start:mat_fin,1)=var_list(:,1);
                        phi_sq_vec(mat_start:mat_fin,1)=var_list.^2;
                        var_matlab_loop=var(var_list);
                        mat_write_stack(i1,:,i2)=[L_in/length_norm,length_norm,avg_a,avg_a_cubed,n_in,var_matlab_loop/norm,...
                            mean(var_list),var_matlab_loop/mean(var_list),skewness(var_list),kurtosis(var_list)];
                    end
                    %Now we compute the 'corrected two pass variance' from Numerical
                    %recipes. It is also called 'compensated variant' on wikipedia.
                    %This is what Matlab uses to calculate the variance
                    var_matlab=var(phi_vec);
                    vec_write(i1,:)=[L_in/length_norm,length_norm,avg_a,avg_a_cubed,n_in,var_matlab/norm,...
                        mean(phi_vec),var_matlab_loop/mean(phi_vec),skewness(phi_vec),kurtosis(phi_vec)];
                    text_path1=strcat(file_path,'\real space averages\variance');
%                     filename=['\variance_circle_windows_avgellnorm ' data_type{dat} '_' type ' ' num2str(width) ' imi_' num '_avg.txt'];
                    filename=['\variance_circle_windows_avgellnorm_interiorOnlypad ' data_type{dat} '_' type ' imi_' num '_avg.txt'];
                    dlmwrite(strcat(text_path1,filename),vec_write(1:i1,:),'Delimiter',',','Precision','%1.15e','newline','pc');
                    i1
                end
                if runs==1
                    figure(1);
                    loglog(2*vec_write(:,1),vec_write(:,6),'Linewidth',2)
                    hold on;
                    loglog(2*vec_write(:,1),avg_a./(pi*(vec_write(:,1)*length_norm).^2),'--r','Linewidth',3)
                else
                    figure(2);
                    loglog(2*vec_write(:,1),vec_write(:,6),'Linewidth',2)
                    hold on;
                    loglog(2*vec_write(:,1),pi*vec_write(:,1).^2,'--r','Linewidth',3)
                end
            end  
%             figure(2);
%             semilogx(2*vec_write(:,1),vec_write(:,8),'Linewidth',2)
%             hold on;
%             figure(3);
%             semilogx(2*vec_write(:,1),vec_write(:,9),'Linewidth',2)
%             hold on;
%             plot([1e-3,1E3],[0,0],'--r','Linewidth',3)
%             figure(4);
%             semilogx(2*vec_write(:,1),vec_write(:,10),'Linewidth',2)
%             hold on;
%             plot([1e-3,1E3],[3,3],'--r','Linewidth',3)
        end
    end
end






