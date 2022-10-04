% write_publication_plots_rescale
%this program will read in data sets and renormalize the x axis

chiq_path='E:\Chieco\Hyperuniformity\HU foams\fourier space data\voronoi data';

type_chi={'Poisson','Einstein_square_delta026','Halton'};
type_mat={'poisson','Einstein_square_delta026','Halton'};
pattern_vec={'pointsINI'};

L_sys_x=1;%sqrt(numel(areas));
L_sys_y=1;%sqrt(numel(areas));
A_sys=L_sys_x*L_sys_y;


for i1=1:numel(pattern_vec)
    pattern=pattern_vec{i1};
    for i2=1:numel(type_mat)
        chiq=dlmread([chiq_path '\Chi_q ' type_chi{i2} ' ' pattern '_points_mono N_cells_ 500000 L_sys_32768.txt']);
        area_read_path='E:\Chieco\Hyperuniformity\HU foams\voronoi\voro_centroids\x_y_area_N lists';
        data=dlmread([area_read_path '\' type_mat{i2} ' centroid_xyAN N_cells_500000 run_1 .txt']);
        bub_area=sum(data(:,3).^2)/sum(data(:,3));
        bub_length_scale=sqrt(bub_area);        
        rho_mat=numel(data(:,3))/A_sys;
        
        chiq_rescale=[chiq(:,1:2)*sqrt(rho_mat)*bub_length_scale,chiq(:,3:4)];
        
        loglog(chiq(2:end,1),chiq(2:end,3),'LineWidth',2);
        hold on
        loglog(chiq_rescale(2:end,1),chiq_rescale(2:end,3),'LineWidth',2);
        
        
        filewrite_name=[chiq_path '\Chi_q_sqrtavganorm ' type_mat{i2} ' ' pattern ' N_cells_500000 L_sys_32768.txt'];
        dlmwrite(filewrite_name,chiq_rescale,'newline','pc')
    end
end
plotedit on