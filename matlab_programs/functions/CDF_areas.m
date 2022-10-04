readPath='D:\Chieco\Hyperuniformity\HU foams\track data';
im_read=[125,250,500];

%we read in the data and show it collapses for A/<A>
figure(1)
for i1=1:numel(im_read)
    bub_lengths_mat=[-1,-1,-1];
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
    phi_avg_a=sum(bub_keep(:,6).^2)/sum(bub_keep(:,6));
    mean_a=mean(bub_keep(:,6));
    norm_a=mean_a;
    i1
    %     [phi_avg_a,mean(bub_keep(:,6))]
    %     [std(bub_keep(:,6))/mean(bub_keep(:,6)),(std(bub_keep(:,6))/mean(bub_keep(:,6)))/sqrt(numel(bub_keep(:,6)))]
    %     [mean(bub_keep(:,1)),std(bub_keep(:,1)),std(bub_keep(:,1))/sqrt(numel(bub_keep(:,1)))]
    numel(bub_keep(:,6))
    a_sq_mean=mean(bub_keep(:,6).^2);
    a_mean_sq=mean(bub_keep(:,6))^2;
    dist_percent_err=sqrt((std(bub_keep(:,6).^2)/a_sq_mean)^2+(2*mean(bub_keep(:,6))*std(bub_keep(:,6))/a_mean_sq)^2)/sqrt(numel(bub_keep(:,6)));
    [a_sq_mean/a_mean_sq,(a_sq_mean/a_mean_sq)*dist_percent_err]
    [mean(bub_keep(:,1)),mean(bub_keep(:,1))/sqrt(numel(bub_keep(:,6))),std(bub_keep(:,1)),std(bub_keep(:,1))/sqrt(numel(bub_keep(:,6)))]
    histogram(bub_keep(:,6)/norm_a);
    if i1==1
        area_vec=bub_keep(:,6)/norm_a;
        area_vec_reg=bub_keep(:,6);
        num_vec=[bub_keep(:,1)];
        std_root_a_vec=std(bub_keep(:,6))/mean(bub_keep(:,6));
    else
        area_vec=[area_vec;bub_keep(:,6)/norm_a];
        area_vec_reg=[area_vec_reg;bub_keep(:,6)];
        num_vec=[num_vec;bub_keep(:,1)];
        std_root_a_vec=[std_root_a_vec;std(bub_keep(:,6))/mean(bub_keep(:,6))];
    end    
    hold on
    cdf_imi_name=['D:\Chieco\Hyperuniformity\HU foams\distribution data\foam data\foam_area_CDF imi_' im_num_str '.txt'];
    cdf_imi=[0,0;sort(bub_keep(:,6)/norm_a),(1:numel(bub_keep(:,6)))'];
    dlmwrite(cdf_imi_name,[cdf_imi(:,1),1-cdf_imi(:,2)/numel(cdf_imi(:,2))],'Newline','pc')
end
numel(area_vec)
a_sq_mean=mean(area_vec_reg.^2);
a_mean_sq=mean(area_vec_reg)^2;
dist_percent_err=sqrt((std(area_vec_reg.^2)/a_sq_mean)^2+(2*mean(area_vec_reg)*std(area_vec_reg)/a_mean_sq)^2)/numel(area_vec_reg)^(1/2);
[a_sq_mean/a_mean_sq,(a_sq_mean/a_mean_sq)*dist_percent_err]
[mean(num_vec),mean(num_vec)/sqrt(numel(num_vec)),std(num_vec),std(num_vec)/sqrt(numel(num_vec))]

cdf_mat=[0,0;sort(area_vec),(1:numel(area_vec))'];
histogram(area_vec)
output=histogram(cdf_mat(2:end,1),50,'Normalization','cdf');
cdf_ys=output.Values;
cdf_xs=output.BinEdges;


comp_exp_fit=@(coeff,xs_in)(exp(-(gamma(1+1/coeff(1)).*xs_in).^coeff(1)));
poiss_cdf_fit=@(coeff,xs_in)(coeff(1)*exp(-(gamma(1+1/coeff(1)).*xs_in).^coeff(1)));
a_fit=1;
xs_fit=cdf_xs(1:end-1);
ys_fit=1-cdf_ys;
y_weights=1./abs(ys_fit(ys_fit>0).^2);
alpha=nlinfit(xs_fit(ys_fit>0)',ys_fit(ys_fit>0)',comp_exp_fit,a_fit,'Weights',y_weights(y_weights>0)')
alpha=1.37

figure(2)
semilogy(cdf_mat(:,1),1-cdf_mat(:,2)/numel(cdf_mat(:,2)),'b','LineWidth',2);
hold on
xs=(0:0.01:10);
ys=exp(-(gamma(1+1/alpha(1))*xs).^alpha(1));
plot(xs,ys,'--r','LineWidth',2)
plot(xs,exp(-xs),':k','LineWidth',2)

total_name='D:\Chieco\Hyperuniformity\HU foams\distribution data\foam data\foam_area_CDF total.txt';
dlmwrite(total_name,[cdf_mat(:,1),1-cdf_mat(:,2)/numel(cdf_mat(:,2))],'Newline','pc')

fit_name='D:\Chieco\Hyperuniformity\HU foams\distribution data\foam data\foam_area_CDF comp_exp_fit.txt';
dlmwrite(fit_name,[xs',ys'],'Newline','pc')

exp_name='D:\Chieco\Hyperuniformity\HU foams\distribution data\foam data\foam_area_CDF exponential.txt';
dlmwrite(exp_name,[xs',exp(-xs')],'Newline','pc')

%new we perform the same analysis but for the voronoi patterns.
readPath='D:\Chieco\Hyperuniformity\HU foams\voronoi';
pattern={'Halton','Einstein_square_delta026','poisson','Einstein_square'};

keyboard
for i1=1:3%numel(pattern)
    centroid_dist_mat=[-1,-1,-1];
    points_dist_mat=[-1,-1,-1];
    pat=pattern{i1};
    cell_mat=dlmread([readPath '\voro_centroids\all info list\' pat ' ID_NA_xyp_xyc_vertID N_cells_500000 run_1 .txt']);
    vert_mat=dlmread([readPath '\voro_vertices\' pat ' vertex_xyID N_cells_500000 run_1.txt']);
    %we need to remove the bubbles that connect to the edge. We find the
    %edge vertices and find which bubbles they are on
    vert_left=vert_mat(vert_mat(:,1)==0,3);
    vert_right=vert_mat(vert_mat(:,1)==1,3);
    vert_bot=vert_mat(vert_mat(:,2)==0,3);
    vert_top=vert_mat(vert_mat(:,2)==1,3);
    bad_verts=unique([vert_left;vert_right;vert_top;vert_bot]);
    %now we want to find the bubbles that have a bad vertex
    bub_ids=cell_mat(:,1);
    vert_ids_mat=cell_mat(:,8:end);
    sz_of_verts=size(vert_ids_mat);
    for i2=1:numel(bad_verts)
        bub_id_list=find(vert_ids_mat==bad_verts(i2));
        [bub_rows,bub_cols]=ind2sub(sz_of_verts,bub_id_list);
        bub_ids(bub_rows)=-1;
    end
    bub_keep=cell_mat(bub_ids>0,:);
    area_vec=bub_keep(:,3)/mean(bub_keep(:,3));
%     figure(1)
%     histogram(area_vec);
    pat    
    numel(bub_keep(:,1));
    stat_num=sqrt(numel(bub_keep(:,3)));
    a_sq_mean=mean(bub_keep(:,3).^2);
    a_mean_sq=mean(bub_keep(:,3))^2;
    dist_percent_err=sqrt((std(bub_keep(:,3).^2)/a_sq_mean)^2+(2*mean(bub_keep(:,3))*std(bub_keep(:,3))/a_mean_sq)^2)/sqrt(numel(bub_keep(:,3)));
    [a_sq_mean/a_mean_sq,(a_sq_mean/a_mean_sq)*dist_percent_err]
    [mean(bub_keep(:,2)),mean(bub_keep(:,2))/sqrt(numel(bub_keep(:,2))),std(bub_keep(:,2)),std(bub_keep(:,2))/sqrt(numel(bub_keep(:,2)))]
    
% %     [std(bub_keep(:,3))/mean(bub_keep(:,3)),std(bub_keep(:,3))/mean(bub_keep(:,3))/sqrt(numel(bub_keep(:,3)))]
%     [mean(bub_keep(:,2)), mean(bub_keep(:,2))/stat_num,std(bub_keep(:,2)), std(bub_keep(:,2))/stat_num]
    hold on
    figure(2)
    cdf_mat=[0,0;sort(area_vec),(1:numel(area_vec))'];
    semilogy(cdf_mat(:,1),1-cdf_mat(:,2)/numel(cdf_mat(:,2)),'LineWidth',2);
    
    %This lets us fit a compressed exponential to the data but I do not
    %think this is the appropriate comparrison
%     figure(3)
%     histogram(area_vec)
%     output=histogram(cdf_mat(2:end,1),150,'Normalization','cdf');
%     cdf_ys=output.Values;
%     cdf_xs=output.BinEdges;    
%     close Figure 3
%     
%     a_fit=2;
%     xs_fit=cdf_xs(1:end-1);
%     ys_fit=1-cdf_ys;
%     y_weights=1./abs(ys_fit(ys_fit>0));
%     alpha=nlinfit(xs_fit(ys_fit>0)',ys_fit(ys_fit>0)',comp_exp_fit,a_fit,'Weights',y_weights(y_weights>0)')
%     
%     xs=(0:0.01:4);
%     ys=exp(-(gamma(1+1/alpha(1))*xs).^alpha(1));
%     plot(xs,ys,'--r','LineWidth',2)
    
end
   