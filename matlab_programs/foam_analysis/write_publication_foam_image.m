nums={'125','250','500'};
mult=[2,1,0];

%This is for a full size image
x_corner=128;
y_corner=2779;

lin_scale=22; %the number of pixels in 1 mm
x_o=x_corner+108;
y_o=y_corner;
y_f=y_o-60;

im_size=2728;

% % %This is for a zoom in image
% x_corner=1261;
% y_corner=2304;
% 
% lin_scale=22; %the number of pixels in 1 mm
% x_o=x_corner+53;
% y_o=y_corner-38;
% y_f=y_o-22;
% 
% im_size=1100;


weak_fit=[500,0;250,0;125,0];
for n_count=1:numel(nums)
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
    
    box_size=13;
    %Now we construct our circular mask
    x1 = box_size/2;
    y1 = box_size/2;
    radius = box_size/2;
    
    [xx,yy] = meshgrid((1:box_size)-y1,(1:box_size)-x1);
    mask=zeros(box_size,box_size)+1.0;
    mask((xx.^2 + yy.^2)<radius^2)=0.0;
    mask=double(mask);
    
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
            x_f=floor(x_o+5*d_scale); 
%             sqrt(avga)/(21.8)
%             im_size=10*d_scale;
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
        
        figure('pos',[20 50 1120 1120])
        ax1 = axes('Position',[0.025 0.025 0.95 0.95]);
        axis equal
        axis off
        hold on
        fig_imi(y_f:y_o,x_o:x_f)=0;        
%         imshow(fig_imi(y_corner-im_size:y_corner,x_corner:x_corner+im_size)/255); 
        places_plot=places(and(and(places(:,1)>=x_corner,places(:,1)<=x_corner+im_size),...
                                and(places(:,2)>=y_corner-im_size,places(:,2)<=y_corner)),:);
        plot(places_plot(:,1)-x_corner,places_plot(:,2)-y_corner+im_size,'o','Color','k',...
                  'LineWidth',2,'MarkerSize',7);
              
        scalebar=[x_o-x_corner,abs((y_o+y_f)/2-y_corner+im_size-max(places_plot(:,2)));
                  x_f-x_corner,abs((y_o+y_f)/2-y_corner+im_size-max(places_plot(:,2)))];
        plot(scalebar(:,1),scalebar(:,2),'k','LineWidth',12)        
        
        fig = gcf;
        filesave=['E:\Chieco\My Papers\hyperuniformity foam\images\unedited images\foam image plots\foam_points_' tag ' imi_' num '.emf'];
   
        
%         fig = gcf;
%         filesave=['E:\Chieco\My Papers\hyperuniformity foam\images\unedited images\foam image plots\foam_' tag '_scale_abar_zoom imi_' num '.emf'];

        exportgraphics(fig,filesave,'ContentType','vector')
        close Figure 1  
    end
end

