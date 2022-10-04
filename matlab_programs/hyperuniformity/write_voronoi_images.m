
centroid_path='F:\Chieco\Hyperuniformity\HU anneal\particle positions\voro_centroids';
pattern='poisson_neighbor_anneal';
cell_txt='3000';
rescale=1;

im_path='F:\Chieco\Hyperuniformity\HU anneal\images';

for i1=0:36
    centroid_info=[centroid_path '\x_y_area_N lists\' pattern ' cenneighbors_xyAN N_cells_' cell_txt ' step_' num2str(i1) '.txt'];    
    my_try=dlmread(centroid_info);
    
%     figure(3)
%     h_areas=my_try(:,3);
%     histogram(h_areas/mean(h_areas))
%     hold on
%     sum(areas.^2)/sum(areas)
%     numel(areas)
%     mean(areas)
%     figure(4)
%     hold on
%     histogram(my_try(:,4))
    
    [v_pos,connections]=voronoin(my_try(:,1:2));
    [vx,vy]=voronoi(my_try(:,1),my_try(:,2));    
    
    figure('pos',[20 50 920 920])
    ax1 = axes('Position',[0.025 0.025 0.95 0.95]);
    hold on
    plot(my_try(:,1),my_try(:,2),'or');
%     plot(v_pos(:,1),v_pos(:,2),'ok')
    axis equal
    xlim([0 rescale])
    ylim([0 rescale])
    plot(vx,vy,'b')

    
    fig = gcf;    
    filesave=[im_path '\' pattern ' voronoi_centers_images step_' num2str(i1) '.png'];
    
    fig.PaperPositionMode
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 8];
    print(fig,'-r400','-dpng',filesave);
    close Figure 1     
    
end