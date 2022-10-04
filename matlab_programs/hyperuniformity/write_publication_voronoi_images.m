%new we perform the same analysis but for the voronoi patterns.
readPath='E:\Chieco\Hyperuniformity\HU foams\voronoi';
pattern={'poisson','Einstein_square_delta026','Halton'};

imi_path='E:\Chieco\My Papers\hyperuniformity foam\images\unedited images';

color_mat=[204/255,0,0,;0,102/255,0;0,0,204/255];

for i1=1:numel(pattern)
    centroid_dist_mat=[-1,-1,-1];
    points_dist_mat=[-1,-1,-1];
    pat=pattern{i1};
    cell_mat=dlmread([readPath '\voro_centroids\all info list\' pat ' ID_NA_xyp_xyc_vertID N_cells_500000 run_1 .txt']);
%     vert_mat=dlmread([readPath '\voro_vertices\x_y_area lists\' pat ' vert_xyA N_cells_5000 run_1 .txt']);
%     arc_mat=dlmread([readPath '\voro_films\x_y_area lists\' pat ' arc_xyA N_cells_5000 run_1 .txt']);
    for i2=1:2
        if i2==1
            voro_mat=cell_mat(:,4:5);
            plot_mat=cell_mat(:,4:5);  
            area_vec=cell_mat(:,3);
            dist_scale=sqrt(sum(area_vec.^2)/sum(area_vec));
            eqsize_scale=5*dist_scale;
            type='points';
        elseif i2==2
            plot_mat=cell_mat(:,6:7);
            area_vec=cell_mat(:,3);
            dist_scale=sqrt(sum(area_vec.^2)/sum(area_vec));
            type='centroids';
        elseif i2==3
            plot_mat=vert_mat(:,1:2);
            dist_scale=sqrt(1/numel(vert_mat(:,1)));
            type='verts';
        else
            plot_mat=arc_mat(:,1:2);
            area_vec=arc_mat(:,3);
            dist_scale=sum(area_vec.^2)/sum(area_vec);
            type='arcs';
        end
        plot_keep_x=plot_mat(and(plot_mat(:,1)>=0.5-(eqsize_scale+dist_scale),...
                                 plot_mat(:,1)<=0.5+eqsize_scale+dist_scale),:);
        clear plot_mat
        plot_keep=plot_keep_x(and(plot_keep_x(:,2)>=0.5-(eqsize_scale+dist_scale),...
                                  plot_keep_x(:,2)<=0.5+eqsize_scale+dist_scale),:);
        clear plot_keep_x
                
        if i2==1
            figure('pos',[20 50 920 920])
            ax1 = axes('Position',[0.025 0.025 0.95 0.95]);
            hold on            
            [v_pos,connections]=voronoin(plot_keep(:,1:2));
            [vx,vy]=voronoi(plot_keep(:,1),plot_keep(:,2));
            plot(plot_keep(:,1),plot_keep(:,2),'o','Color',color_mat(i1,:),...
                 'MarkerFaceColor',color_mat(i1,:),'LineWidth',1,'MarkerSize',6);
            axis equal
            xlim([0.5-eqsize_scale 0.5+eqsize_scale])
            ylim([0.5-eqsize_scale 0.5+eqsize_scale])
            plot(vx,vy,'Color',color_mat(i1,:),'LineWidth',2)
        else
            plot(plot_keep(:,1),plot_keep(:,2),'o','Color',color_mat(i1,:),'LineWidth',2,'MarkerSize',7);
        end
    end
    scalebar=[0.5-eqsize_scale+dist_scale/2,0.5-eqsize_scale+dist_scale/2;
              0.5-eqsize_scale+3*dist_scale/2,0.5-eqsize_scale+dist_scale/2];
              
    plot(scalebar(:,1),scalebar(:,2),'k','LineWidth',12) 
    
    fig = gcf;
    set(gca,'visible','off')
    filesave=[imi_path '\' pat ' voronoi_' type '_images.emf'];
    exportgraphics(fig,filesave,'ContentType','vector')
%     fig.PaperPositionMode
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 8 8];
%     print(fig,'-r600','-dpng',filesave);
    close Figure 1
end
   