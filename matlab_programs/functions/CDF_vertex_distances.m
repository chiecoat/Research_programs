readPath='E:\Chieco\Hyperuniformity\HU foams\patterns_to_analyze';
im_read=[125,250,500];

%we read in the data and show it collapses for A/<A>
figure(1)
for i1=1:numel(im_read)
    bub_lengths_mat=[-1,-1,-1];
    im_num_str=num2str(im_read(i1));  
    %the arc lengths are the same as the distances between vertices. 
    stats_mat=dlmread([readPath '\arc_stats_mat imi_' im_num_str '.txt']);
    %we need to remove vertces that are along the edge
    vert_left=stats_mat(stats_mat(:,1)>min(stats_mat(:,1)),:);
    vert_right=vert_left(vert_left(:,1)<max(vert_left(:,1)),:);
    vert_bot=vert_right(vert_right(:,2)>min(stats_mat(:,2)),:);
    vert_in=vert_bot(vert_bot(:,2)<max(vert_bot(:,2)),:);    
    dists_keep=vert_in(:,3);
    histogram(dists_keep/mean(dists_keep));
    i1
    %     numel(dists_keep)
    %     [std(dists_keep)/mean(dists_keep),std(dists_keep)/mean(dists_keep)/sqrt(numel(dists_keep))]
    numel(dists_keep)
    dist_sq_mean=mean(dists_keep.^2);
    dist_mean_sq=mean(dists_keep)^2;
    dist_percent_err=sqrt((std(dists_keep.^2)/dist_sq_mean)^2+(2*mean(dists_keep)*std(dists_keep)/dist_mean_sq)^2)/sqrt(numel(dists_keep));
    [dist_sq_mean/dist_mean_sq,(dist_sq_mean/dist_mean_sq)*dist_percent_err]
    hold on
    if i1==1
        length_vec=dists_keep/mean(dists_keep);
        length_reg=dists_keep;
    else
        length_vec=[length_vec;dists_keep/mean(dists_keep)];
        length_reg=[dists_keep;length_reg];
    end    
    cdf_imi_name=['E:\Chieco\Hyperuniformity\HU foams\distribution data\foam data\foam_vert_distance_CDF imi_' im_num_str '.txt'];
    cdf_imi=[0,0;sort(dists_keep/mean(dists_keep)),(1:numel(dists_keep))'];
    dlmwrite(cdf_imi_name,[cdf_imi(:,1),1-cdf_imi(:,2)/numel(cdf_imi(:,2))],'Newline','pc')
    cdf_actual_name=['E:\Chieco\Hyperuniformity\HU foams\distribution data\foam data\foam_vert_distance_CDF_Actual imi_' im_num_str '.txt'];
    dlmwrite(cdf_actual_name,[cdf_imi(2:end,1),cdf_imi(2:end,2)/numel(cdf_imi(:,2))],'Newline','pc')
end

numel(length_reg)
dist_sq_mean=mean(length_reg.^2);
dist_mean_sq=mean(length_reg)^2;
dist_percent_err=sqrt((std(length_reg.^2)/dist_sq_mean)^2+(2*mean(length_reg)*std(length_reg)/dist_mean_sq)^2)/sqrt(numel(length_reg));
[dist_sq_mean/dist_mean_sq,(dist_sq_mean/dist_mean_sq)*dist_percent_err]

cdf_mat=[0,0;sort(length_vec),(1:numel(length_vec))'];
numel(length_vec)
[std(length_vec)/mean(length_vec),std(length_vec)/mean(length_vec)/sqrt(numel(length_vec))]
histogram(length_vec)
output=histogram(cdf_mat(2:end,1),150,'Normalization','cdf');
cdf_ys=output.Values;
cdf_xs=output.BinEdges;


comp_exp_fit=@(coeff,xs_in)(exp(-(gamma(1+1/coeff(1)).*xs_in).^coeff(1)));
a_fit=1.5;
xs_fit=cdf_xs(1:end-1);
ys_fit=1-cdf_ys;
y_weights=1./abs(ys_fit(ys_fit>0));
alpha=nlinfit(xs_fit(ys_fit>0)',ys_fit(ys_fit>0)',comp_exp_fit,a_fit,'Weights',y_weights(y_weights>0)');

figure(2)
% semilogy(cdf_mat(:,1),1-cdf_mat(:,2)/numel(cdf_mat(:,2)),'b','LineWidth',2);
loglog(cdf_mat(:,1),cdf_mat(:,2)/numel(cdf_mat(:,2)),'b','LineWidth',2);
hold on
xs=(0:0.01:2);
ys=exp(-(gamma(1+1/alpha(1))*xs).^alpha(1));
plot(xs,ys,'--r','LineWidth',2)
plot(xs,exp(-xs),':k','LineWidth',2)

total_name='E:\Chieco\Hyperuniformity\HU foams\distribution data\foam data\foam_vert_distance_CDF total.txt';
dlmwrite(total_name,[cdf_mat(:,1),1-cdf_mat(:,2)/numel(cdf_mat(:,2))],'Newline','pc')

total_name='E:\Chieco\Hyperuniformity\HU foams\distribution data\foam data\foam_vert_distance_CDF_actual total.txt';
dlmwrite(total_name,[cdf_mat(2:end,1),cdf_mat(2:end,2)/numel(cdf_mat(:,2))],'Newline','pc')

%new we perform the same analysis but for the voronoi patterns.
readPath='E:\Chieco\Hyperuniformity\HU foams\voronoi\voro_films\x_y_area lists';
pattern={'Halton','Einstein_square_delta026','poisson','Einstein_square'};


for i1=1:3%numel(pattern)
    centroid_dist_mat=[-1,-1,-1];
    points_dist_mat=[-1,-1,-1];
    pat=pattern{i1};
    stats_mat=dlmread([readPath '\' pat ' arc_xyA N_cells_500000 run_1.txt']);
    %we need to remove the vertices that lie only along the edge
    vert_left=stats_mat(stats_mat(:,1)>0,:);
    vert_right=vert_left(vert_left(:,1)<1,:);
    vert_bot=vert_right(vert_right(:,2)>0,:);
    vert_in=vert_bot(vert_bot(:,2)<1,:);
    dists_keep=vert_in(:,3);
    length_vec=dists_keep/mean(dists_keep);
    figure(1)
    histogram(length_vec);
    pat
    numel(length_vec)
    dist_sq_mean=mean(length_vec.^2);
    dist_mean_sq=mean(length_vec)^2;
    dist_percent_err=sqrt((std(length_vec.^2)/dist_sq_mean)^2+(2*mean(length_vec)*std(length_reg)/dist_mean_sq)^2)/sqrt(numel(length_vec));
    [dist_sq_mean/dist_mean_sq,(dist_sq_mean/dist_mean_sq)*dist_percent_err]
    %     numel(length_vec)
%     [std(dists_keep)/mean(dists_keep),std(dists_keep)/mean(dists_keep)/sqrt(numel(length_vec))]
    hold on
    figure(2)
    cdf_mat=[0,0;sort(length_vec),(1:numel(length_vec))'];
%     semilogy(cdf_mat(:,1),1-cdf_mat(:,2)/numel(cdf_mat(:,2)),'LineWidth',2);
    loglog(cdf_mat(:,1),cdf_mat(:,2)/numel(cdf_mat(:,2)),'LineWidth',2);
end
