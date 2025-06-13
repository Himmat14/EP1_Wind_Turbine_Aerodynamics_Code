clc
clear

blade_data = readmatrix("forces_blade.xlsx");

R = blade_data(:,1);
c = blade_data(:,2);
twist = blade_data(:,3);
cl = blade_data(:,8);
cd = blade_data(:,9);
aoa = blade_data(:,5);
a = blade_data(:,10);
a_prime = blade_data(:,12);


theta = twist + aoa;

dR=R(2)-R(1);

V = 10;

rho = 1.225;

omega = (3000/60);

W = sqrt(((1-a).^2).*V^2+((R.*omega).^2).*(1+a_prime.^2));

dF_z = rho.*(W.^2).*(cl.*cos(theta)+cd.*sin(theta));

dF_t = rho.*(W.^2).*(cl.*sin(theta)-cd.*cos(theta));

F_z = sum(dF_z.*dR);

F_T = sum(dF_t).*dR;


