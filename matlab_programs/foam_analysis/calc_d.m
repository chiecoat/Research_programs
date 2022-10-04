%This is a shor tprogram just to we don't clutter the foam track program.
%We will calculate the value "d" for the square 3D printed trough, given a
%fill height for the foaming solution. The fill height must be in mm!
function d_height=calc_d(h_fill)

%First we need to know the parameters of the trough, it is a square that
%is 336mmx336mm 
h_total=10.1; %mm
h_gap=2.4; %mm
h_trough=7.7; %mm
w_trough=10; %mm


L_long_trough=336; %mm
L_short_trough=316; %mm
L_avg_trough=326; %mm;

%there is an outlet that can be used to remove any unwanted object in the
%trough. ITs size is 40.5mm x 10mm
outlet_width=10; %mm
outlet_length=40.5; %mm


%The volume of the fluid that fills the container
if h_fill<=w_trough
    v_fluid=h_fill*L_long_trough*h_total;
else
    v_fluid=h_fill*L_long_trough*h_total-(h_fill-w_trough)*L_short_trough*h_trough;
end

%The height of the fluid in the container
h_fluid=v_fluid/(4*w_trough*L_avg_trough+outlet_width*outlet_length);

d_height=(h_trough+h_gap/2)-h_fluid;

end
