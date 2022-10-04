% +
% NAME: Chi_q_calc_sparse_image
%
% PURPOSE:
%     This program will find the spectral density for bidisperse mixture of
%     disks. To do so it will call the "partial structure" code. Patrial
%     structure can be found in the "Function" folder in MatLab Programs.
%
% CATEGORY:
%     Hyperuniformity
%
% CALLING SEQUENCE:
%    FFT_my_image
%
% INPUTS: 
%    There are no inputs but we need to read in particle positoins for
%    whatever configuration we are trying to analyze.
%
% OPTIONAL INPUTS: (none)
%
% KEYWORD PARAMETERS: (none)
%
% OUTPUTS: Chi_q
%    The function will write data files that are the spectral density 
%    versus wave number for our input configuration.     
%
% 
% SIDE EFFECTS: (none)
%
% MODIFICATION HISTORY:
%    written by: Remi Dreyfus, UPenn
%     edited by: A. Chieco, UPenn 2017
%-
% clear all 
% close all

read_path='F:\Chieco\Hyperuniformity\HU foams\voronoi images';
savePath='F:\Chieco\Hyperuniformity\HU foams';
imiPathvec={'voro_centroids','voro_films','voro_vertices'};
structure={'centroids','films','vertex'};
imi_nums=[2644,1382,1009];
n_runs=10;
rescale=2848;

f_path=[savePath '\fourier space data'];
for i1=1:numel(imiPathvec)
    imi_folder=imiPathvec{i1};
    feature=structure{i1};
    for i2=1:numel(imi_nums)
        n_bubs=imi_nums(i2);
        for i3=1:n_runs
            imiPathTot=[read_path '\' imi_folder '\voronoi_' feature ' nbubs=' num2str(imi_nums(i2)) ' rescale=' num2str(rescale) ' imi_' num2str(i3) '.png'];
            feature_imi=imread(imiPathTot);
            
            feature_imi(feature_imi<150)=0;
            feature_imi(feature_imi>150)=1;
            
            L_sys_x=numel(feature_imi(1,:));
            L_sys_y=numel(feature_imi(:,1));
            A_sys=L_sys_x*L_sys_y;
            
            areas=feature_imi(feature_imi>0);
            
            avg_a=sum(areas.^2)/(A_sys);
            avg_a_cubed=sum(areas.^4)/(A_sys);
            N_tot=numel(feature_imi(feature_imi>0));
            
            a_nums=unique(feature_imi);
            
            imi_norm=1/sqrt(sum(a_nums.^2));
            imi_in=feature_imi;
            clear feature_imi
            spec_add=fft_radial_large(imi_in,imi_norm);
            clear imi_in
            spec_tot=[spec_add(:,1)*sqrt(avg_a)/sqrt(A_sys),spec_add(:,2)];
            %     loglog(spec_tot(2:end,1),spec_tot(2:end,2))%'Color', [255, 215, 0]./255)
            %     hold on
            bidi_name=['fftshift_binary Structure_f_voronoi_' feature ' nbubs=' num2str(imi_nums(i2)) ' rescale=' num2str(rescale) ' imi_' num2str(i3) '.txt'];
            dlmwrite([f_path '\' bidi_name],spec_tot,'Precision',12,'newline','pc');
            bark=0;
        end
    end
end