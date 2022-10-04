%This code works for bidisperse packings of disks. It has not been
%generalized for arbitrary polydispersity but could be.
function [n_phi] = circle_window_overlapexact_periodic(sys_corners,xy_tot,num,side,point_pattern)

%We find the radius of our circle
rad=side;

r_small=min(xy_tot(:,3));
d_small_inside=sqrt(rad^2-r_small^2);
r_large=max(xy_tot(:,3));
d_large_inside=sqrt(rad^2-r_large^2);

mat_small=xy_tot(xy_tot(:,3)==r_small,:);
mat_large=xy_tot(xy_tot(:,3)==r_large,:);

%we set the dimensions of our matrix first
xdim=sys_corners(2,1)-sys_corners(1,1);
ydim=sys_corners(2,2)-sys_corners(1,2);

%We find the limits for x
xmin=sys_corners(1,1);
xmax=sys_corners(2,1);
%We find the limits for y
ymin=sys_corners(1,2);
ymax=sys_corners(2,2);
c_points=[xmin,ymin;xmin,ymax;xmax,ymin;xmax,ymax];

%we find the corners of our windows from a random number generator
if or(strcmp(point_pattern,'Low_disc'),strcmp(point_pattern,'Lattice'))==0
    win=rand(num,2);
    winx=xmin+(xmax-xmin)*win(:,1);
    winy=ymin+(ymax-ymin)*win(:,2);
end

%n_phi will hold the number of centers identified when we want to
%calculate the standard deviation.
n_phi=zeros(1,num);
for i1=1:num
    %now we find how many centers we have in the data
    cen=[winx(i1),winy(i1)];
    %c_points are the locations of the corner of our system. If our circle
    %lies within sqrt(2)*rad of the corners then there is some overlap for
    %the circular window and the system boundary
    c_ov=c_points(sqrt((c_points(:,1)-cen(1)).^2+(c_points(:,2)-cen(2)).^2)<=sqrt(2)*rad,:);
    if numel(c_ov)>0
        %we check to see which corner we are closest to
        c_test=c_points(sqrt((c_points(:,1)-cen(1)).^2+(c_points(:,2)-cen(2)).^2)==...
                        min(sqrt((c_points(:,1)-cen(1)).^2+(c_points(:,2)-cen(2)).^2)),:);         
        cen_prime=[xdim*(c_points(:,1)-c_test(1)),ydim*(c_points(:,2)-c_test(2))];
        cen_array=[cen(1)+cen_prime(:,1),cen(2)+cen_prime(:,2)];
    else
        %Now we check if our circle does not overlap any of the edges. If
        %that is the case then the window exists in the bulk of the system and
        %we can continue
        no_ov=and(and(cen(1)-rad >= xmin,cen(1)+rad < xmax), ...
            and(cen(2)-rad >= ymin,cen(2)+rad < ymax));
        if no_ov==1
            cen_array=cen;
        else  
            %we check to see if our circular window overlaps only the minimum x boundary
            if and(cen(1)-rad < xmin,and(cen(2)-rad >= ymin,cen(2)+rad < ymax))==1
                cen_array=[cen;cen(1)+xdim,cen(2)];
            end
            %we check to see if our circular window overlaps only the maximum x boundary
            if and(cen(1)+rad >= xmax,and(cen(2)-rad >= ymin,cen(2)+rad < ymax))==1
                cen_array=[cen;cen(1)-xdim,cen(2)];
            end
            %we check to see if our circular window overlaps only the minimum y boundary
            if and(cen(2)-rad < ymin,and(cen(1)-rad >= xmin,cen(1)+rad < xmax))==1
                cen_array=[cen;cen(1),cen(2)+ydim];
            end
            %we check to see if our circular window overlaps only the maximum y boundary
            if and(cen(2)+rad >= ymax,and(cen(1)-rad >= xmin,cen(1)+rad < xmax))==1
                cen_array=[cen;cen(1),cen(2)-ydim];
            end
        end
    end
    sum_vec=zeros(3,numel(cen_array(:,1)));
    for i3=1:numel(cen_array(:,1))
        %we find the distance between the center of our observation window
        %and the circles that have any part that lies within the circle
        %If the window size is smaller than the circle size then the entire
        %window may lie inside of a particle. If the window size is larger
        %than the particle then we have to check for entire particles within 
        %the window.
        if rad<r_small
            dist_small_tot=sqrt((mat_small(:,1)-cen_array(i3,1)).^2+(mat_small(:,2)-cen_array(i3,2)).^2);
            dist_small=dist_small_tot(dist_small_tot+rad>r_small);
            %This identifies windows that lie entirely within a small
            %particle
            area_total_small=mat_small(dist_small_tot+rad<=r_small,1:3);
            dist_small_entire=dist_small_tot(dist_small_tot+rad<=r_small);            
            dist_large_tot=sqrt((mat_large(:,1)-cen_array(i3,1)).^2+(mat_large(:,2)-cen_array(i3,2)).^2);
            dist_large=dist_large_tot(dist_large_tot+rad>r_large);
            %This identifies windows that lie entirely within a small
            %particle
            area_total_large=mat_large(dist_large_tot+rad<=r_large,1:3);
            dist_large_entire=dist_large_tot(dist_large_tot+rad<=r_large);
            %We also need to know if the entire window lies within a circle,
            %then the are covered is the area of the entire circle
            sum_vec(1,i3)=numel(dist_small_entire)*pi*rad^2;
            sum_vec(2,i3)=numel(dist_large_entire)*pi*rad^2;
        end
        %This is if the window is larger than the small particles but
        %smaller than the large particles
        if and(rad<r_large,rad>r_small)==1
            dist_small=sqrt((mat_small(:,1)-cen_array(i3,1)).^2+(mat_small(:,2)-cen_array(i3,2)).^2);
            dist_large_tot=sqrt((mat_large(:,1)-cen_array(i3,1)).^2+(mat_large(:,2)-cen_array(i3,2)).^2);
            dist_large=dist_large_tot(dist_large_tot+rad>r_large);
            %Here we find the small particles that lie entirely within the
            %window
            area_total_small=mat_small(dist_small<=rad-r_small,1:3);
            sum_vec(1,i3)=numel(area_total_small(:,1))*pi*r_small^2;
            %This identifies windows that lie entirely within a small
            %particle
            area_total_large=mat_small(dist_large_tot+rad<=r_large,1:3);
            dist_large_entire=dist_large_tot(dist_large_tot+rad<=r_large);
            %We also need to know if the entire window lies within a circle,
            %then the are covered is the area of the entire circle
            sum_vec(2,i3)=numel(dist_large_entire)*pi*rad^2;
        end
        if rad>r_large
            dist_small=sqrt((mat_small(:,1)-cen_array(i3,1)).^2+(mat_small(:,2)-cen_array(i3,2)).^2);
            dist_large=sqrt((mat_large(:,1)-cen_array(i3,1)).^2+(mat_large(:,2)-cen_array(i3,2)).^2);
            %Here we find the small particles that lie entirely within the
            %window
            area_total_small=mat_small(dist_small<=rad-r_small,1:3);
            sum_vec(1,i3)=numel(area_total_small(:,1))*pi*r_small^2;
            %Here we find the large particles that lie entirely within the
            %window
            area_total_large=mat_large(dist_large<=rad-r_large,1:3);
            sum_vec(2,i3)=numel(area_total_large(:,1))*pi*r_large^2;     
        end 
        %There are three types of particles that overlap with our circular
        %window whose radius is rad.
        %If the any part of the particle lies within out circular
        %observation window then we can determine the overlapping area
        area_in_small=mat_small(and(dist_small>rad-r_small,dist_small<rad+r_small),1:3);
        dist_in_small=dist_small(and(dist_small>rad-r_small,dist_small<rad+r_small));
        areas_overlap_small=r_small^2*acos((r_small^2+dist_in_small.^2-rad^2)./(2*r_small*dist_in_small))+...
                            rad^2*acos((rad^2+dist_in_small.^2-r_small^2)./(2*rad*dist_in_small))-...
                            0.5*sqrt((-dist_in_small+r_small+rad).*(dist_in_small+r_small-rad).*(dist_in_small-r_small+rad).*(dist_in_small+r_small+rad));
        area_in_large=mat_large(and(dist_large>rad-r_large,dist_large<rad+r_large),1:3);  
        dist_in_large=dist_large(and(dist_large>rad-r_large,dist_large<rad+r_large));
        areas_overlap_large=r_large^2*acos((r_large^2+dist_in_large.^2-rad^2)./(2*r_large*dist_in_large))+...
                            rad^2*acos((rad^2+dist_in_large.^2-r_large^2)./(2*rad*dist_in_large))-...
                            0.5*sqrt((-dist_in_large+r_large+rad).*(dist_in_large+r_large-rad).*(dist_in_large-r_large+rad).*(dist_in_large+r_large+rad));
        sum_vec(3,i3)=sum(areas_overlap_small)+sum(areas_overlap_large);
        %This is only to visualize if we have identified the correct
        %circles.
%         circles(cen_array(i3,1),cen_array(i3,2),rad,'FaceColor','none','EdgeColor','k','LineWidth',2);
%         hold on
%         axis equal
%         if numel(area_total_small)>0
%             circles(area_total_small(:,1),area_total_small(:,2),area_total_small(:,3),'FaceColor','b','EdgeColor','b','LineWidth',0.001);
%         end
%         if numel(area_total_large)>0
%             circles(area_total_large(:,1),area_total_large(:,2),area_total_large(:,3),'FaceColor','r','EdgeColor','r','LineWidth',0.001);
%         end
%         if numel(area_out_small)>0
%             circles(area_out_small(:,1),area_out_small(:,2),area_out_small(:,3),'FaceColor','g','EdgeColor','g','LineWidth',0.001);
%         end
%         if numel(area_out_large)>0
%             circles(area_out_large(:,1),area_out_large(:,2),area_out_large(:,3),'FaceColor',[1 0.333 0],'EdgeColor',[1 0.333 0],'LineWidth',0.001);
%         end
%         if numel(area_in_small(:,1))>0
%             circles(area_in_small(:,1),area_in_small(:,2),area_in_small(:,3),'FaceColor','c','EdgeColor','c','LineWidth',0.001);
%         end
%         if numel(area_in_large)>0
%             circles(area_in_large(:,1),area_in_large(:,2),area_in_large(:,3),'FaceColor','m','EdgeColor','m','LineWidth',0.001);
%         end
%         circles(cen_array(i3,1),cen_array(i3,2),rad,'FaceColor','none','EdgeColor','k','LineWidth',2);
%         bark=0;
    end   
    bark=0;
    if numel(sum_vec(1,:))>1
       n_phi(i1)=double(sum(sum(sum_vec)))/double(pi*rad^2);
    else 
        n_phi(i1)=double(sum(sum_vec))/double(pi*rad^2);
    end   
end
n_phi=n_phi';

end