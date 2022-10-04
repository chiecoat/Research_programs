% +
% NAME: voro_plus_plus_file_make
%
% PURPOSE:
%     This program will create a file that can be read by voro++. These
%     files need to be in formatted at id x y z (r), the id must be an 
%     integer and the r is optional. 
%
% CATEGORY:
%     Hyperuniformity
%
% CALLING SEQUENCE:
%    find_h_length
%
% INPUTS: None at the moment as it is a stand alone program, but we need to
% know the filenames.
%
% OPTIONAL INPUTS: (none)
%
% KEYWORD PARAMETERS: (none)
%
% OUTPUTS: Will write a text file for voro++

% SIDE EFFECTS: (none)
%
% MODIFICATION HISTORY:
%    written by: A. Chieco, UPenn, December 2016
%
%-

read_path='D:\Chieco\Hyperuniformity\HU Jammed Packings from Ning Xu';
quench_type='\quenchN1E6';
tol='\tol1.0e-16';

phis=[850];
n_loops=1;
%This is the rati betweens large and small diameters of our disks
x_min=-0.5; x_max=0.5; side_x=x_max-x_min;
y_min=-0.5; y_max=0.5; side_y=y_max-y_min;
A_sys=1; N_tot=1E6; sig_con=1.4;
N_big=N_tot/2; N_small=N_tot/2;

write_path='D:\Chieco\Hyperuniformity\voronoi_constructions\quenchVaryRate\N=1E6\particle_positions_for_voro';
%If we have periodic boundary conditions for our particle positions we want
%to layer a few particles on either side of the boundary. This way we can
%get 
periodic=1;

for i1=1:numel(phis)
    phi_txt=num2str(phis(i1));
    p_frac=phis(i1);
    phi=p_frac/10^(ceil(log10(p_frac+1)));
    %Now we compute radius of a small particle, sig_s
    sig_s=sqrt(8*phi*A_sys/(N_tot*pi*(1+sig_con^2)))/2;
    sig_b=sig_con*sig_s;
    a_small=pi*(sig_s)^2;
    a_big=pi*(sig_b)^2;
    %here is the the phi weighted average area of a particle, <a>
    phi_big=N_big*a_big/A_sys; phi_small=N_small*a_small/A_sys;
    avg_a=(phi_big*a_big+phi_small*a_small)/phi;
    for i2=1:n_loops
        pos_path=['\particle positions\phi.' phi_txt '_' num2str(i2)];
        pos_mat=dlmread([read_path quench_type tol pos_path],'\t',[1 0 N_tot 1]);  
        pos_mat=[pos_mat,[zeros(N_tot/2,1)+sig_s;zeros(N_tot/2,1)+sig_b],zeros(N_tot,1)+1];
        if periodic==1 
            pos_mat_long=zeros(9*N_tot,4);
            pos_mat_long(1:N_tot,:)=pos_mat;
            %we stamp the whole system around our current particles. Then
            %we will take all the particles that exist in a larger box so 
            %can enforce our peridoic boundary.
            add_mat=[side_x,0,3;-side_x,0,-3;0,side_y,5;0,-side_y,-5;...
                     side_x,side_y,4;side_x,-side_y,2;-side_x,side_y,-2;-side_x,-side_y,-4];
             for j1=1:numel(add_mat(:,1))
                 pos_mat_long(j1*N_tot+1:(j1+1)*N_tot,:)=[pos_mat(:,1)+add_mat(j1,1),pos_mat(:,2)+add_mat(j1,2),pos_mat(:,3),zeros(N_tot,1)+add_mat(j1,3)];
             end
             %Now we cut away most of the particles we dont want
             pos_x_keep=pos_mat_long(and(pos_mat_long(:,1)<x_max+10*sig_b,pos_mat_long(:,1)>=x_min-10*sig_b),:);
             pos_tot=pos_x_keep(and(pos_x_keep(:,2)<y_max+10*sig_b,pos_x_keep(:,2)>=y_min-10*sig_b),:);
             n_ids=numel(pos_tot(:,1));
             %We want to write two text files out for peridoic boundary
             %conditions. The first is for voro++ because it needs specific
             %formating. The second will have the same format as voro++
             %file only it will also include a tag for where particles
             %landed from applying the peridoic condtions.
             voro_pp_mat=[(1:1:n_ids)',pos_tot(:,1:2),zeros(n_ids,1),pos_tot(:,3:4)];
             file_write=['\voro_phi' phi_txt '_' num2str(i2) '.txt'];
             dlmwrite([write_path file_write],voro_pp_mat(:,1:end-1),'Delimiter','\t','Precision',8)
             file_write=['\voro_phi' phi_txt '_' num2str(i2) '_tag.txt'];
             dlmwrite([write_path file_write],voro_pp_mat,'Delimiter','\t','Precision',8)
        else
            %Without periodic boundary conditions we just needed to tag the
            %particles with *integer* IDs, add a column for z positions
            %(all zeros) and add a column for particle radii
            pos_tot=pos_mat;
            n_ids=numel(pos_tot(:,1));
            voro_pp_mat=[(1:1:n_ids)',pos_tot(:,1:2),zeros(n_ids,1),pos_tot(:,3)];
            file_write=['\voro_phi' phi_txt '_' num2str(i2) '.txt'];
            dlmwrite([write_path file_write],voro_pp_mat,'Delimiter','\t','Precision',8)
        end
    end
end






