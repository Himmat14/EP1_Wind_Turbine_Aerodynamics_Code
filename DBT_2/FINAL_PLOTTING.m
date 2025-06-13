%plotting 
clear
clc

format long

set(groot,'defaultLineLineWidth',2)  %sets graph line width as 2
set(groot,'defaultAxesFontSize',18)  %sets graph axes font size as 18
set(groot,'defaulttextfontsize',18)  %sets graph text font size as 18
set(groot,'defaultLineMarkerSize',8) %sets line marker size as 8
set(groot,'defaultAxesXGrid','on')   %sets X axis grid on 
set(groot,'defaultAxesYGrid','on')   %sets Y axis grid on
set(groot, 'DefaultAxesBox', 'on')   %sets Axes boxes on

data_fin = readmatrix("testdatafinal_plotting.xlsx");

aerodat = readmatrix("Aerofoildat.xlsx");

power_dat = readmatrix("POWERvsTSR.xlsx");


TSR = power_dat(:,1);
test1 = power_dat(:,2);
test2 = power_dat(:,3);
test3 = power_dat(:,4);
test4 = power_dat(:,5);

%  test 1

final_power = figure;
set(final_power,"WindowState","maximized");

hold on 
plot(TSR,test1)
yline(power_dat(1,6),LineWidth=2,Color='r', LineStyle='-.')
title('Power vs TSR')
xlabel('TSR')
ylabel('Power (W)')
legend('Final Design','Betz limit',Location='east')
axis([0 8 0 80])
hold off

saveas(final_power,'Power vs TSR.png');

%chord plot

Chord_dist = figure;
set(Chord_dist,'WindowState','maximized');

R = data_fin(:,1);
hold on
plot(R,data_fin(:,2));
title('Chord Distribution')
xlabel('radial position (m)')
ylabel('Chord length (m)')
axis([0.05 0.25 0.03 0.11])
hold off

saveas(Chord_dist,'Chord Distribution.png');

%twist plot

Twist_dist = figure;
set(Twist_dist,'WindowState','maximized');

hold on
plot(R,data_fin(:,4));
title('Twist Distribution')
xlabel('radial position (m)')
ylabel('Twist angle (°)')
hold off

saveas(Twist_dist,'Twist Distribution.png')

%aoa plot

AOA_des_dist = figure;
set(AOA_des_dist,'WindowState','maximized');

hold on 
plot(R,data_fin(:,6))
title('alpha_d_e_s_i_g_n distribution')
xlabel('radial position (m)')
ylabel('Angle of Attack')
hold off

saveas(AOA_des_dist,'AOA distribution.png')

%cl des plot

CL_des_dist = figure;
set(CL_des_dist,'WindowState','maximized');

hold on 
plot(R,data_fin(:,7))
title('Cl_d_e_s_i_g_n distribution')
xlabel('radial position (m)')
hold off

saveas(CL_des_dist,'CL_design dist.png')

%cl plot

CL_dist = figure;
set(CL_dist,'WindowState','maximized');

hold on 
plot(R,data_fin(:,8))
title('C_L distrubution')
xlabel('radial position (m)')
ylabel('C_L')
hold off

saveas(CL_dist,'CL distribution.png')

%cd plot

CD_dist = figure;
set(CD_dist,'WindowState','maximized');

hold on
plot(R,data_fin(:,9))
title('C_D distribution')
xlabel('radial position (m)')
ylabel('C_D')
hold off

saveas(CD_dist,'CD distribution.png')

%a plot

a_dist = figure;
set(a_dist,'WindowState','maximized');

hold on
plot(R,data_fin(:,10))
yline(0.33)
title('Axial Induction Factor distribution')
xlabel('radial position (m)')
ylabel('a')
hold off

saveas(a_dist,'Axial induction factor.png')

% a_prime plot

a_prime_dist = figure;
set(a_prime_dist,'WindowState','maximized');

hold on
plot(R,data_fin(:,12))
title('Angular Induction Factor distribution')
xlabel('radial position (m)')
ylabel('a prime')
hold off

saveas(a_prime_dist,'Angular induction factor.png')









