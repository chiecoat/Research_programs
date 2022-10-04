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
    for i2=1:1
        if i2==1
            voro_mat=cell_mat(:,4:5);
            plot_mat=cell_mat(:,4:5);  
            area_vec=cell_mat(:,3);
            dist_scale=sqrt(sum(area_vec.^2)/sum(area_vec));
            eqsize_scale=20*dist_scale;
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
            plot(plot_keep(:,1),plot_keep(:,2),'o','Color',color_mat(i1,:),...
                 'MarkerFaceColor',color_mat(i1,:),'LineWidth',1,'MarkerSize',6);
            axis equal
            xlim([0.5-eqsize_scale 0.5+eqsize_scale])
            ylim([0.5-eqsize_scale 0.5+eqsize_scale])
            numel(plot_keep(:,1))
        else
            plot(plot_keep(:,1),plot_keep(:,2),'o','Color',color_mat(i1,:),'LineWidth',2,'MarkerSize',7);
        end
    end
    scalebar=[0.5-eqsize_scale+3*dist_scale/2,0.5-eqsize_scale+3*dist_scale/2;
              0.5-eqsize_scale+13*dist_scale/2,0.5-eqsize_scale+3*dist_scale/2];
              
    plot(scalebar(:,1),scalebar(:,2),'k','LineWidth',12) 
    
    fig = gcf;
    set(gca,'visible','off')
    filesave=[imi_path '\' pat ' points_' type '_images.emf'];
    exportgraphics(fig,filesave,'ContentType','vector')
    close Figure 1
end
% scale_bar_loc=2848*dist_scale;
% % 
% % nums={'125','250','500'};
% % 
% % weak_fit=[500,0;250,0;125,0];
for n_count=1:3%numel(nums)
    num=nums{n_count};    
    
    savePath='E:\Chieco\Hyperuniformity\HU foams\patterns_to_analyze';
    vPath=[savePath '\vert_stats_mat imi_' num];
    vertex_write=[vPath '.txt'];
    vertex_total=dlmread(vertex_write);
    
    L_sys_x=max(vertex_total(:,1))-min(vertex_total(:,1));
    L_sys_y=max(vertex_total(:,2))-min(vertex_total(:,2));
    A_sys=L_sys_x*L_sys_y;
    
    bubPath=[savePath '\bub_stats_mat imi_' num];
    bubbles_write=[bubPath '.txt'];
    bub_total=dlmread(bubbles_write);
    size(bub_total)
    
    arcPath=[savePath '\arc_stats_mat imi_' num];
    arcs_write=[arcPath '.txt'];
    arcs_total=dlmread(arcs_write);
    
    raw_imi=imread(['E:\Chieco\My Papers\hyperuniformity foam\images\unedited images\foam image plots\imi read\raw imi_' num '.png']);
    raw_imi=double(raw_imi);
    arc_imi=imread(['E:\Chieco\My Papers\hyperuniformity foam\images\unedited images\foam image plots\imi read\arcs_overlay imi_' num '.png']);
    arc_imi=double(arc_imi);
    arc_locs=find(arc_imi<10);
    
    for mats=2:2
        if mats==1
            foam_imi=raw_imi;
            fig_imi=raw_imi;
            places=vertex_total(:,1:2);
            fin_value=5;
            tag='vertex';
            L_sys_x=max(places(:,1))-min(places(:,1));
            L_sys_y=max(places(:,2))-min(places(:,2));
            A_sys=L_sys_x*L_sys_y;
            d_scale=sqrt(A_sys/numel(places(:,1)));
            x_f=floor(x_o+d_scale);            
        elseif mats==2
            foam_imi=raw_imi;
            fig_imi=raw_imi;
            places=bub_total(:,1:2);
            fin_value=10;
            tag='bubbles';
            avga=sum(bub_total(:,3).^2)/sum(bub_total(:,3));
            d_scale=sqrt(mean(avga));%sqrt(avga)
            [avga,mean(bub_total(:,3))]/(21.8^2);
            weak_fit(n_count,2)=mean(bub_total(:,3));
            x_f=floor(x_o+d_scale); 
%             sqrt(avga)/(21.8)
            im_size=22.5*d_scale;
        else
            foam_imi=arc_imi;
            fig_imi=raw_imi;
            foam_imi(arc_locs)=1;
            places=arcs_total(:,1:2);
            fin_value=15;
            tag='arcs';
            avg_ell=sum(arcs_total(:,3).^2)/sum(arcs_total(:,3));
            d_scale=avg_ell;
            x_f=floor(x_o+d_scale); 
        end    
        
        plot_keep_x=places(and(places(:,1)>=L_sys_x/2-L_sys_y/2,...
                               places(:,1)<=L_sys_x/2+L_sys_y/2),:);
        clear plot_mat
        plot_keep=plot_keep_x(and(plot_keep_x(:,2)>=0,...
                               plot_keep_x(:,2)<=L_sys_y),:);
        
        figure('pos',[20 50 1120 1120])
        ax1 = axes('Position',[0.025 0.025 0.95 0.95]); 
        hold on     
        plot(plot_keep(:,1),plot_keep(:,2),'o','Color','k',...
                 'MarkerFaceColor','k','LineWidth',1.5,'MarkerSize',6);
        axis equal
        axis off
        
        scalebar=[L_sys_x/2-L_sys_y/2+d_scale,d_scale;
                  L_sys_x/2-L_sys_y/2+6*d_scale,d_scale];
              
        plot(scalebar(:,1),scalebar(:,2),'k','LineWidth',12) 
        
        
        
        fig = gcf;
        filesave=['E:\Chieco\My Papers\hyperuniformity foam\images\unedited images\foam image plots\foam_' tag '_scale_abar_zoom imi_' num '.emf'];

        exportgraphics(fig,filesave,'ContentType','vector')       
        close Figure 1 
        
    end
    
end

