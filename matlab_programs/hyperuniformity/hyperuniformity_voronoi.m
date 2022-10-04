% clear all
% close all
% format long g

file_path='F:\Chieco\Hyperuniformity\HU foams';
pattern='Einstein_square_delta026';
% data_type={[pattern ' vert_xyA'],[pattern ' arc_xyA'],[pattern ' centroid_xyAN'],[pattern ' pointsINI_xyAN']};
% data_folder={'\voronoi\voro_vertices\x_y_area lists','\voronoi\voro_films\x_y_area lists',...
%              '\voronoi\voro_centroids\x_y_area_N lists','\voronoi\voro_centroids\x_y_area_N lists'};
data_type={[pattern ' pointsINI_xyAN'],[pattern ' centroid_xyAN'],[pattern ' vert_xyA'],};
data_folder={'\voronoi\voro_centroids\x_y_area_N lists','\voronoi\voro_centroids\x_y_area_N lists',...
             '\voronoi\voro_vertices\x_y_area lists'};
       
r_start=1; r_fin=1;

%the normalization length for all of the data is just the sqrt(<a>) of the
%bubble.

num_cells_str={'5000'};

for num_cells=1:numel(num_cells_str)    
    for dat=1:1%numel(data_type)
        for n_run=1:r_fin
            n_cells=num_cells_str{num_cells};
            foam_norm=dlmread(['F:\Chieco\Hyperuniformity\HU foams\voronoi'...
                '\voro_centroids\x_y_area_N lists\' pattern ' centroid_xyAN N_cells_' n_cells ' run_' num2str(n_run) ' .txt']);
            area_norm=sum(foam_norm(:,3).^2)/sum(foam_norm(:,3));
            bub_length_norm=sqrt(area_norm);
            f_path=[file_path '\fourier space data\voronoi data'];
            %first we read in each file
            fname=[file_path  data_folder{dat} '\' data_type{dat} ' N_cells_' n_cells ' run_' num2str(n_run) ' .txt'];
            mat=dlmread(fname);
            x_min=0; x_max=1;
            y_min=0; y_max=1;
            L_sys=1;
            A_sys=L_sys^2;
            sys_specs=[x_min,y_min;x_max,y_max];
            mat_keep_tot=mat;
            for runs=1:numel(mat(1,3:end))+1
                if runs==1
                    continue
%                     areas=mat(:,3);
%                     type='points_weighted';
                elseif runs==numel(mat(1,3:end))+1
                    areas=mat(:,3)*0+1;
                    type='points_mono';
                else
                    continue
%                     areas=mat(:,end);
%                     type='centers_coord';
                end
                a_vec=areas;
                phi=sum(areas)/A_sys;
                norm=1;
                avg_a=sum(areas.^2)/sum(areas);
                avg_a_cubed=sum(areas.^4)/sum(areas);
                if runs==1
                    length_norm=sqrt(avg_a);
                else
                    length_norm=sqrt(A_sys/numel(areas));
                end
                %we find how many windows and hat values here
                c_rads1=logspace(log10(0.5E-2*length_norm),log10(length_norm/2),15);
                c_rads2=logspace(log10(length_norm/2),log10(0.3*L_sys),30);
                c_rads=[c_rads1,c_rads2(2:end)];
                N_windows_vec=zeros(1,numel(c_rads))+1E5;%round(L_sys^2./(c_rads.^2));
                %Now we test for hyperuniformity at a given box size L and keep all of
                %the raw data in a matrix
                i1_start=1; i1_fin=numel(c_rads);
                vec_write=zeros(i1_fin,10);
                mat_write_stack=zeros(i1_fin,10);
                for i1=1:i1_fin
                    n_in=N_windows_vec(i1);
                    L_in=c_rads(i1);
                    n_phi_vec=r_fin*n_in;
                    phi_vec=zeros(n_phi_vec,1);
                    phi_sq_vec=phi_vec;
                    for i2=r_start:r_fin
                        var_list=circle_window_points(sys_specs,[mat_keep_tot(:,1:2)],n_in,L_in);
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
                    text_path1=strcat(file_path,'\real space averages\variance\voronoi data');
                    filename=['\variance_circle_windows ' data_type{dat} '_' type ' voro_nbubs_' n_cells '_run_' num2str(n_run) '.txt'];
                    dlmwrite(strcat(text_path1,filename),vec_write(1:i1,:),'Delimiter',',','Precision','%1.15e','newline','pc');
                    i1
                end
                figure(1)                
                loglog(2*vec_write(:,1),vec_write(:,6),'Linewidth',2)      
                hold on
                bark=0;
                %         text_path2=strcat(file_path,'\real space data\loop data');
                %         for i3=r_start:r_fin
                %             var_mat_write=mat_write_stack(:,:,i3);
                %             filename=['\variance ' data_type{dat} '_' type ' voro_nbubs_' nbubs_str '_loop_',num2str(i3),'.txt'];
                %             dlmwrite(strcat(text_path2,filename),var_mat_write,'Delimiter',',','Precision','%1.15e','newline','pc');
                %         end
            end
        end
        loglog(2*vec_write(:,1),(pi*(vec_write(:,1)).^2),'--r','Linewidth',3)        
        figure(2);
        semilogx(2*vec_write(:,1),vec_write(:,8),'Linewidth',2)
        hold on; 
        figure(3); 
        semilogx(2*vec_write(:,1),vec_write(:,9),'Linewidth',2)
        hold on; 
        plot([1e-3,1E3],[0,0],'--r','Linewidth',3)
        figure(4); 
        semilogx(2*vec_write(:,1),vec_write(:,10),'Linewidth',2)
        hold on; 
        plot([1e-3,1E3],[3,3],'--r','Linewidth',3)
        keyboard
    end
end
    
    
    
    
    
