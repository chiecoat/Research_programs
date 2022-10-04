% +
% NAME: elongation_for_hyper_foam
%
% PURPOSE: 
%   This program will read in different patterns and find the elongation of
%   each of the cells in that pattern
%
% CATEGORY:
%   data analysis
%
% CALLING SEQUENCE:
%   none this will be a standalone program
%
% INPUTS: 
%   non but we will need to put in the filenames and folders
%
% OPTIONAL INPUTS: (none)
%
% KEYWORD PARAMETERS: (none)
%
% OUTPUTS: Will write files to destination folders
% 
% SIDE EFFECTS: (none)
%
% MODIFICATION HISTORY:
%    written by: A. Chieco UPenn, September 2020
%-


readPath='E:\Chieco\Hyperuniformity\HU foams\track data';
im_read=[125,250,500];

%we read in the data to calculate the elongations for each foam images
max_vec_check=[0,0,0];
for i1=1:numel(im_read)  
    im_num_str=num2str(im_read(i1));  
    vert_mat=dlmread([readPath '\vertices\vertices imi_' im_num_str '_update.txt']);
    bub_mat=dlmread([readPath '\bubbles\bubbles imi_' im_num_str '_final.txt']);  
    stats_mat=dlmread([readPath '\stats\stats_bub_nxy_pace imi_' im_num_str '_final.txt']);
    %we need to remove the bubbles that connect to the edge. We do this by
    %identifying vertices that touch cell=1 and remove any bubbles with
    %that vertex from the list
    cell_ids=vert_mat(:,end-2:end);
    id_list=find(cell_ids==1);
    [vert_rows,vert_cols]=ind2sub(size(cell_ids),id_list);
    bub_ids=bub_mat(:,1);
    vert_ids_mat=bub_mat(:,5:end);
    sz_of_verts=size(vert_ids_mat);
    for i2=1:numel(vert_rows)
        bub_id_list=find(vert_ids_mat==vert_rows(i2));
        [bub_rows,bub_cols]=ind2sub(sz_of_verts,bub_id_list);
        bub_ids(bub_rows)=-1;        
    end
    bub_keep=stats_mat(bub_ids>0,:);
    bub_elong=bub_keep(:,10);
    i1
    E_mean=mean(bub_elong);
    E_mean_sq=E_mean^2;
    E_sq_mean=mean(bub_elong.^2);    
%   E_sq_mean,E_mean_sq,var(bub_elong),mean(bub_elong)]
%   dist_percent_err=sqrt((std(bub_elong.^2)/E_sq_mean)^2+(2*mean(bub_elong)*std(bub_elong)/E_mean_sq)^2)/sqrt(numel(bub_elong));
%   [E_sq_mean/E_mean_sq,(E_sq_mean/E_mean_sq)*dist_percent_err]
    E_mean_err=std(bub_elong)/sqrt(numel(bub_elong));
    E_sq_mean_err=std(bub_elong.^2)/sqrt(numel(bub_elong));
    [E_mean,10*E_mean_err,E_sq_mean,10*E_sq_mean_err]
%     
    figure(1)
    output=histogram(bub_elong,50,'Normalization','probability');
    edges=output.BinEdges;
    plot_xs=(edges(1:end-1)+edges(2:end))/2;
    plot_ys=output.Values;
    close Figure 1
    figure(3)
    plot(plot_xs,plot_ys,'LineWidth',3)
    hold on
    max_vec_check(i1)=max(bub_elong);
    if i1==1
        elong_vec=bub_elong;
        full_foam_mat=[bub_keep(:,2:3),bub_keep(:,1),bub_keep(:,10),bub_keep(:,12)*0+i1];
    else
        elong_vec=[bub_elong;elong_vec];
        full_foam_mat=[full_foam_mat;bub_keep(:,2:3),bub_keep(:,1),bub_keep(:,10),bub_keep(:,12)*0+i1];
    end    
    figure(2)
    cdf_imi=[0.5,0;sort(bub_elong),(1:numel(bub_elong))'];
    plot(cdf_imi(:,1),1-cdf_imi(:,2)/numel(cdf_imi(:,2)))
    hold on
%     cdf_imi_name=['E:\Chieco\Hyperuniformity\HU foams\distribution data\foam data\foam_elong_CDF imi_' im_num_str '.txt'];
%     dlmwrite(cdf_imi_name,[cdf_imi(:,1),1-cdf_imi(:,2)/numel(cdf_imi(:,2))],'Newline','pc')
end
full_foam_mat=full_foam_mat(full_foam_mat(:,4)<min(max_vec_check),:);
elong_vec=elong_vec(elong_vec<min(max_vec_check));

'foam all'
E_mean=mean(elong_vec);
E_mean_sq=E_mean^2;
E_sq_mean=mean(elong_vec.^2);
%   E_sq_mean,E_mean_sq,var(bub_elong),mean(bub_elong)]
%   dist_percent_err=sqrt((std(bub_elong.^2)/E_sq_mean)^2+(2*mean(bub_elong)*std(bub_elong)/E_mean_sq)^2)/sqrt(numel(bub_elong));
%   [E_sq_mean/E_mean_sq,(E_sq_mean/E_mean_sq)*dist_percent_err]
E_mean_err=std(elong_vec)/sqrt(numel(elong_vec));
E_sq_mean_err=std(elong_vec.^2)/sqrt(numel(elong_vec));
[E_mean,10*E_mean_err,E_sq_mean,10*E_sq_mean_err]

% histogram(elong_vec)
figure(1)
output=histogram(elong_vec,50,'Normalization','probability');
edges=output.BinEdges;
plot_xs=(edges(1:end-1)+edges(2:end))/2;
plot_ys=output.Values;
close Figure 1
figure(3)
plot(plot_xs,plot_ys,'LineWidth',3)
plot([mean(elong_vec),mean(elong_vec)],[0,10])

figure(4)
cdf_full=[0.5,0;sort(elong_vec),(1:numel(elong_vec))'];
plot(cdf_full(:,1),1-cdf_full(:,2)/numel(cdf_full(:,2)))
hold on

[elong_vals,ids]=sort(full_foam_mat(:,4));
elong_tail=flip(ids);
foam_look=full_foam_mat(elong_tail,:);

[full_foam_mat_sort,sort_ids]=sort(full_foam_mat(:,4));
full_foam_mat=[full_foam_mat(sort_ids,:),(1:numel(elong_vec))'];
color_mat=[204,0,0; 255,140,0; 0,153,0; 0,204,204; 0,0,204; 139,0,255]/255;
n_sides_vec=[4,5,6,7,8,9];


num_vec_max=logspace(0,log10(1.8),4);
num_vec_min=logspace(log10(1/num_vec_max(3)),0,3);
num_scale_vec=[num_vec_min,num_vec_max(2:4)];
for plot_col=1:6
    if plot_col<6
        cdf_val_foam=full_foam_mat(full_foam_mat(:,3)==n_sides_vec(plot_col),6);
        elong_n_foam=full_foam_mat(full_foam_mat(:,3)==n_sides_vec(plot_col),4);
    else
        cdf_val_foam=[1;full_foam_mat(full_foam_mat(:,3)>=n_sides_vec(plot_col),6)];
        elong_n_foam=[0;full_foam_mat(full_foam_mat(:,3)>=n_sides_vec(plot_col),4)];
    end    
    figure(4)
    y_scale=num_scale_vec(plot_col);
    loglog(elong_n_foam,(1-cdf_val_foam/numel(full_foam_mat(:,2)))*y_scale,'.','Color',color_mat(plot_col,:))
    hold on
end
xlim([1 3.5]);
ylim([5E-5 5]);
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.03],'LineWidth',3,'FontSize',10)


%new we perform the same analysis but for the voronoi patterns.
readPath='E:\Chieco\Hyperuniformity\HU foams\voronoi';
pattern={'Halton','Einstein_square_delta026','poisson','Einstein_square'};
for i2=1:3    
    % for quik and dirty group meting slides, to be deleted later
    elong_mat=dlmread(['E:\Chieco\Hyperuniformity\HU foams\voronoi\' pattern{i2} ' ID_NA_perim_elong N_cells_500000 run_1.txt']);
    elong_vec=elong_mat(:,5);
    pattern{i2}
    E_mean=mean(elong_vec);
    E_mean_sq=E_mean^2;
    E_sq_mean=mean(elong_vec.^2);
    %   E_sq_mean,E_mean_sq,var(bub_elong),mean(bub_elong)]
    %   dist_percent_err=sqrt((std(bub_elong.^2)/E_sq_mean)^2+(2*mean(bub_elong)*std(bub_elong)/E_mean_sq)^2)/sqrt(numel(bub_elong));
    %   [E_sq_mean/E_mean_sq,(E_sq_mean/E_mean_sq)*dist_percent_err]
    E_mean_err=std(elong_vec)/sqrt(numel(elong_vec));
    E_sq_mean_err=std(elong_vec.^2)/sqrt(numel(elong_vec));
    [E_mean,10*E_mean_err,E_sq_mean,10*E_sq_mean_err]
    figure(1)
    output=histogram(elong_vec,100,'Normalization','probability');
    edges=output.BinEdges;
    plot_xs=(edges(1:end-1)+edges(2:end))/2;
    plot_ys=output.Values;
    close Figure 1
    figure(3)
    plot(plot_xs,plot_ys,'LineWidth',3)    
    hold on
    plot([mean(elong_vec),mean(elong_vec)],[0,10])
    figure(4)
    cdf_full=[0.5,0;sort(elong_vec),(1:numel(elong_vec))'];
    plot(cdf_full(:,1),1-cdf_full(:,2)/numel(cdf_full(:,2)))
    
    [full_voro_mat_sort,sort_ids]=sort(elong_mat(:,5));
    full_voro_mat=[elong_mat(sort_ids,:),(1:numel(sort_ids))'];
    
    num_vec_max=logspace(0,log10(1.8),4);
    num_vec_min=logspace(log10(1/num_vec_max(3)),0,3);
    num_scale_vec=[num_vec_min,num_vec_max(2:4)];
    for plot_col=1:6
        if plot_col<6
            cdf_val_voro=full_voro_mat(full_voro_mat(:,2)==n_sides_vec(plot_col),6);
            elong_n_voro=full_voro_mat(full_voro_mat(:,2)==n_sides_vec(plot_col),5);
        else
            cdf_val_voro=[1;full_voro_mat(full_voro_mat(:,2)>=n_sides_vec(plot_col),6)];
            elong_n_voro=[0;full_voro_mat(full_voro_mat(:,2)>=n_sides_vec(plot_col),5)];
        end
        figure(4)
        y_scale=num_scale_vec(plot_col);
        plot(elong_n_voro,(1-cdf_val_voro/numel(elong_mat(:,5)))*y_scale,'.','Color',color_mat(plot_col,:))
        hold on
    end
end
keyboard

read_path='E:\Chieco\My Papers\hyperuniformity foam\images\unedited images\elongation PDF';
set(gca,'YTickLabel',[],'XTickLabel',[])
fig = gcf;
filesave=[read_path '\elongation_foams.emf'];

%  
% for i1=3:-1:1%numel(pattern)
%     centroid_dist_mat=[-1,-1,-1];
%     points_dist_mat=[-1,-1,-1];
%     pat=pattern{i1};
%     cell_mat=dlmread([readPath '\voro_centroids\all info list\' pat ' ID_NA_xyp_xyc_vertID N_cells_500000 run_1 .txt']);
%     vert_mat=dlmread([readPath '\voro_vertices\' pat ' vertex_xyID N_cells_500000 run_1.txt']);
%     %we need to remove the bubbles that connect to the edge. We find the
%     %edge vertices and find which bubbles they are on
%     vert_left=vert_mat(vert_mat(:,1)==0,3);
%     vert_right=vert_mat(vert_mat(:,1)==1,3);
%     vert_bot=vert_mat(vert_mat(:,2)==0,3);
%     vert_top=vert_mat(vert_mat(:,2)==1,3);
%     bad_verts=unique([vert_left;vert_right;vert_top;vert_bot]);
%     %now we want to find the bubbles that have a bad vertex
%     bub_ids=cell_mat(:,1);
%     vert_ids_mat=cell_mat(:,8:end);
%     sz_of_verts=size(vert_ids_mat);
%     for i2=1:numel(bad_verts)
%         bub_id_list=find(vert_ids_mat==bad_verts(i2));
%         [bub_rows,bub_cols]=ind2sub(sz_of_verts,bub_id_list);
%         bub_ids(bub_rows)=-1;
%     end
%     %we have eliminated all of the edge bubbles and want to calculate the
%     %elongations of the remaining ones
%     bub_keep=cell_mat(bub_ids>0,:);
%     area_vec=bub_keep(:,3);
%     perim_vec=zeros(numel(bub_keep(:,3)),1);    
%     elong_vec=zeros(numel(bub_keep(:,3)),1);
%     for i3=1:numel(bub_keep(:,1))
%         n_sides=bub_keep(i3,2);
%         bub_verts=bub_keep(i3,8:8+n_sides-1);
%         perim_sum_vec=zeros(n_sides,1);
%         for i4=1:n_sides
%             if i4<n_sides
%                 v1=bub_verts(i4);
%                 v2=bub_verts(i4+1);
%             else
%                 v1=bub_verts(i4);
%                 v2=bub_verts(1);
%             end
%             vert_point_1=vert_mat(vert_mat(:,3)==v1,1:2);
%             vert_point_2=vert_mat(vert_mat(:,3)==v2,1:2);
%             dist=sqrt(sum((vert_point_1-vert_point_2).^2));
%             perim_sum_vec(i4)=dist;
%         end
%         perim=sum(perim_sum_vec);
%         perim_vec(i3)=perim;
%         elong_vec(i3)=perim/sqrt(4*pi*area_vec(i3));    
%         if mod(i3,1000)==0
%             i3
%             bub_data_mat=[bub_keep(:,1:3),perim_vec,elong_vec];
%         end        
%     end
%     bub_data_mat=[bub_keep(:,1:3),perim_vec,elong_vec];
%     cell_data_name=['E:\Chieco\Hyperuniformity\HU foams\voronoi\' pat ' ID_NA_perim_elong N_cells_500000 run_1.txt'];
%     dlmwrite(cell_data_name,bub_data_mat,'Newline','pc')    
% end
%    
% 
