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

figure(2)
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
%we make a histogram for the PDF
figure(1)
output=histogram(elong_vec,31,'Normalization','pdf');
edges=output.BinEdges;
plot_xs=(edges(1:end-1)+edges(2:end))/2;
plot_ys=output.Values;
plot_xs=plot_xs(plot_ys~=0);
plot_ys=plot_ys(plot_ys~=0);
close Figure 1
figure(4)
loglog(plot_xs,plot_ys,'LineWidth',8,'color','k')
hold on
bin_edge_vec=(edges(1):edges(2)-edges(1):2*edges(end));

for plot_col=1:6
    full_foam_for_n=full_foam_mat;
    if plot_col<6
        full_foam_for_n(full_foam_mat(:,3)~=n_sides_vec(plot_col),4)=max(bin_edge_vec);
        figure(1)
        output=histogram(full_foam_for_n(:,4),bin_edge_vec,'Normalization','pdf');
        edges=output.BinEdges;
        plot_xs=(edges(1:end-1)+edges(2:end))/2;
        plot_ys=output.Values;
        close Figure 1
    else
        full_foam_for_n(full_foam_mat(:,3)<n_sides_vec(plot_col),4)=max(bin_edge_vec);
        figure(1)
        output=histogram(full_foam_for_n(:,4),bin_edge_vec,'Normalization','pdf');
        edges=output.BinEdges;
        plot_xs=(edges(1:end-1)+edges(2:end))/2;
        plot_ys=output.Values;
        close Figure 1
    end    
    plot_xs=plot_xs(plot_ys~=0);
    plot_ys=plot_ys(plot_ys~=0);
    figure(4)
    loglog(plot_xs(1:end-1),plot_ys(1:end-1),'LineWidth',4,'Color',color_mat(plot_col,:))
    hold on
end
xlim([1 3.5]);
% ylim([5E-5 5]);
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.03],'LineWidth',3,'FontSize',10)


%new we perform the same analysis but for the voronoi patterns.
readPath='E:\Chieco\Hyperuniformity\HU foams\voronoi';
pattern={'Halton','Einstein_square_delta026','poisson','Einstein_square'};
bin_nums=[18,30,50];
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
    figure(2)
    cdf_full=[0.5,0;sort(elong_vec),(1:numel(elong_vec))'];
    plot(cdf_full(:,1),1-cdf_full(:,2)/numel(cdf_full(:,2)))
    
    [full_voro_mat_sort,sort_ids]=sort(elong_mat(:,5));
    full_voro_mat=[elong_mat(sort_ids,:),(1:numel(sort_ids))'];
    
    figure(1)
    output=histogram(elong_vec,bin_nums(i2),'Normalization','pdf');
    edges=output.BinEdges;
    plot_xs=(edges(1:end-1)+edges(2:end))/2;
    plot_ys=output.Values;
    plot_xs=plot_xs(plot_ys~=0);
    plot_ys=plot_ys(plot_ys~=0);
    close Figure 1
    figure(4)
    loglog(plot_xs,plot_ys,'LineWidth',8,'color','k')
    hold on
    bin_edge_vec=(edges(1):edges(2)-edges(1):2*edges(end));
    for plot_col=1:6
        full_voro_for_n=full_voro_mat;
        if plot_col<6
            full_voro_for_n(full_voro_mat(:,2)~=n_sides_vec(plot_col),5)=max(bin_edge_vec);
            figure(1)
            output=histogram(full_voro_for_n(:,5),bin_edge_vec,'Normalization','pdf');
            edges=output.BinEdges;
            plot_xs=(edges(1:end-1)+edges(2:end))/2;
            plot_ys=output.Values;
            close Figure 1
        else
            full_voro_for_n(full_voro_mat(:,2)<n_sides_vec(plot_col),5)=max(bin_edge_vec);
            figure(1)
            output=histogram(full_voro_for_n(:,5),bin_edge_vec,'Normalization','pdf');
            edges=output.BinEdges;
            plot_xs=(edges(1:end-1)+edges(2:end))/2;
            plot_ys=output.Values;
            close Figure 1
        end
        plot_xs=plot_xs(plot_ys~=0);
        plot_ys=plot_ys(plot_ys~=0);
        figure(4)
        loglog(plot_xs(1:end-1),plot_ys(1:end-1),'LineWidth',4,'Color',color_mat(plot_col,:))
        hold on
    end
end
keyboard

read_path='E:\Chieco\My Papers\hyperuniformity foam\images\unedited images\elongation PDF';
set(gcf,'pos',[20 50 1000 700*(3.06/4.76)])
set(gca,'XMinorTick','off','YMinorTick','on','TickLength',[0.03, 0.015],'LineWidth',2)
set(gca,'YTickLabel',[],'XTickLabel',[])
fig = gcf;
filesave=[read_path '\elongation_foams.emf'];
exportgraphics(fig,filesave,'ContentType','vector')
