% +
% NAME: dry_foam_moments
%
% PURPOSE: this will dtermine various moments of dry foam distributions      
%
% CATEGORY:
%     Dry foam analysis
%
% CALLING SEQUENCE:
%    dry_foam_moments
%
% INPUTS: 
%    None
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
%    written by: A. Chieco, UPenn, June 2021
%-

format long

% % % %first thing we want to do is read in a data file
% % % read_folder='E:\Chieco\Dry Foam Moments';
% % % %this reads data from different files, they track with dry foam data from
% % % %the 2d tracking folder
% % % date='0711_nomin';
% % % 
% % % %the data rows are time (min), A (mm^2), dAdt (mm^2/min), n_sides, id
% % % dry_foam_data=readmatrix([read_folder '\' date '_dry_foam_data_t_A_dAdt_n.xlsx']);
% % % %this is a relic of the way the data was written, currently the first row
% % % %is all 0s. Can rewrite data and remove this
% % % dry_foam_data=dry_foam_data(2:end,:);
% % % 
% % % times=(min(dry_foam_data(:,1)):max(dry_foam_data(:,1)));
% % % ids=unique(dry_foam_data(:,end));
% % % 
% % % %first we plot time sequence data fr <<n>>_p and <A^p>/<A>^p
% % % n_sides_vec=(4:8);
% % % powers=(1:20);
% % % n_times=numel(times);
% % % n_powers=numel(powers);
% % % A_moment_mat=zeros(n_times,n_powers);
% % % avg_avg_n_data_mat=zeros(n_times,n_powers);
% % % for i1=1:n_times
% % %     bubs_t_i1=dry_foam_data(dry_foam_data(:,1)==i1,:);
% % %     area_tot=sum(bubs_t_i1(:,2));
% % %     avg_a=mean(bubs_t_i1(:,2));
% % %     for i2=1:n_powers
% % %         power=powers(i2);
% % %         avg_a_power=mean(bubs_t_i1(:,2).^power);
% % %         A_moment_mat(i1,i2)=avg_a_power/(avg_a^power);
% % %         %now we need to calculate <<n>>_p
% % %         n_data_vec=zeros(1,numel(n_sides_vec));
% % %         for i3=1:numel(n_sides_vec)
% % %             n_sides=n_sides_vec(i3);
% % %             n_sides_data=bubs_t_i1(bubs_t_i1(:,end-1)==n_sides,:);
% % %             n_data_vec(i3)=n_sides*sum(n_sides_data(:,2))/area_tot;
% % %         end
% % %         avg_avg_n_data_mat(i1,i2)=sum(n_data_vec);   
% % %     end    
% % % end
% % % 
% % % figure(1)
% % % plot(times,A_moment_mat(:,2),'-o')
% % % hold on
% % % figure(2)
% % % plot(times,avg_avg_n_data_mat(:,1),'-o')
% % % hold on

%keyboard

% clear all 
% close all

%there is a problem with the tracking data. We have to remove bubbles from
%face switching so there are many missing bubbles in any frame. We can
%compare the data to data from just raw bubble information that is not
%tracked

% This is the folder used for the PRE
folder_titles='2019-07-11 nitrogen_fill_height=73mm dry perimeter';
num_start=52100;
loop_start=100;
loop_end=1089;
%First we read in all of the files with the bubbles information
file_path=['E:\Chieco\2d tracking\' folder_titles];
stat_path=[file_path '\bubbles'];
file_name='bubbles imi';

%The conversion from pixel to mm is 1mm=22.005pix  and we convert to real
%units in the when reading in the files.
mm_conv=24.8;

%first we plot time sequence data fr <<n>>_p and <A^p>/<A>^p
n_sides_vec=(4:8);
powers=(1:20);
% n_times=numel(times);
n_powers=numel(powers);

%this is the most bubbles we can have before needing to downsize
%N_max=251;
N_max=240;
x_im_size=3884;
x_mid=x_im_size/2;
y_im_size=2068;
y_mid=y_im_size/2;
x_scale=sqrt(N_max*x_im_size/y_im_size);
y_scale=sqrt(N_max*y_im_size/x_im_size);

overall_scale=sqrt(y_im_size/x_im_size);


N_bub_mat=zeros(loop_end-loop_start+1,3);
for rounds=1:2
    A_moment_raw_data=zeros(loop_end-loop_start+1,n_powers);
    avg_avg_n_raw_data=zeros(loop_end-loop_start+1,n_powers);
    
    for i1=loop_start:loop_end
        file=[stat_path '\' file_name '_' num2str(num_start+i1) ' final.txt'];
        %The bubble information we read in is
        %n_sides,x_bub_cen,y_bub_cen,perimeter,del_perimeter,area,del_area,
        %circularity,del_circularity,elongation,del_elongation, "time". The del
        %indicates uncertainty in that measurment The "time" is just how many
        %frames we are from the initial frame
        bubble_info=dlmread(file);
        %In order for track to work we need to have the coordinates of particle
        %centers first
        bub_stack_all=[bubble_info(:,1:3),bubble_info(:,4)./(mm_conv^2),zeros(numel(bubble_info(:,1)),1)+i1];
        bub_stack_all=bub_stack_all(and(bub_stack_all(:,4)>0,bub_stack_all(:,4)<1E4),:);
        if or(rounds==1,numel(bub_stack_all(:,1))<N_max)==1
            x_min=0;
            x_max=max(bubble_info(:,1));
            y_min=0;
            y_max=max(bubble_info(:,2));
        elseif and(rounds==2,numel(bub_stack_all(:,1))>N_max)==1
            avg_a=mean(bubble_info(bubble_info(:,4)<1E4*mm_conv^2,4));
            avg_r=sqrt(avg_a/pi);
            x_min=x_mid-x_scale*avg_r*sqrt(pi/4);
            x_max=x_mid+x_scale*avg_r*sqrt(pi/4);
            y_min=y_mid-y_scale*avg_r*sqrt(pi/4);
            y_max=y_mid+y_scale*avg_r*sqrt(pi/4);
        end
        bub_stack_x=bub_stack_all(and(bub_stack_all(:,1)>=x_min,bub_stack_all(:,1)<=x_max),:);
        bub_stack=bub_stack_x(and(bub_stack_x(:,2)>=y_min,bub_stack_x(:,2)<=y_max),:);        
        areas=bub_stack(:,4);
        n_sides=unique(bub_stack(:,3));
        n_sides_vec=n_sides(and(n_sides>=3,n_sides<=13));
        area_tot=sum(areas);
        avg_a=mean(areas);
        for i2=1:n_powers
            power=powers(i2);
            avg_a_power=mean(areas.^power);
            A_moment_raw_data(i1-loop_start+1,i2)=avg_a_power/(avg_a^power);
            %now we need to calculate <<n>>_p
            n_data_vec=zeros(1,numel(n_sides_vec));
            for i3=1:numel(n_sides_vec)
                n_sides=n_sides_vec(i3);
                n_sides_data=bub_stack(bub_stack(:,3)==n_sides,:);
                n_data_vec(i3)=(n_sides*sum(n_sides_data(:,4).^power))/(sum(areas.^power));
            end
            avg_avg_n_raw_data(i1-loop_start+1,i2)=sum(n_data_vec);
        end
        N_bub_mat(i1-loop_start+1,1)=i1-loop_start+1;
        N_bub_mat(i1-loop_start+1,1+rounds)=numel(areas);
        if and(rounds==2,i1==loop_start)==1
            figure(2)
            foam_im=imread('E:\foam coarsening\73mm fill height circle trough\2019-07-11 dry_run nitrogen_fill_height=73mm\skeleton images\skeleton DSC_52200.png');
            watershed_mat=watershed(foam_im,8);
            rgb = label2rgb(watershed_mat,'parula',[.5 .5 .5]);
            imshow(rgb);
            hold on
            plot(bub_stack_all(:,1),bub_stack_all(:,2),'ok')
            plot(bub_stack(:,1),bub_stack(:,2),'or')
%             keyboard
            close Figure 2
        end
    end
    
    for i3=1:5
        subplot(1,2,1)
        semilogy((loop_start:loop_end)'-loop_start+1,A_moment_raw_data(:,i3),'-o')
        hold on
        plot([1,loop_end],[mean(A_moment_raw_data(:,i3)),mean(A_moment_raw_data(:,i3))],'--k','Linewidth',2)
        subplot(1,2,2)
        plot((loop_start:loop_end)'-loop_start+1,avg_avg_n_raw_data(:,i3),'-o')
        hold on
        plot([1,loop_end],[mean(avg_avg_n_raw_data(:,i3)),mean(avg_avg_n_raw_data(:,i3))],'--k','Linewidth',2)
    end
    
    
end

