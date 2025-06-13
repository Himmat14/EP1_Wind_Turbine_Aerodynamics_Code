clear all
clc
format long

data_test = readmatrix("FINALDISTRIBUTIONS.csv");

%test 1: NREL S823 , SG6043    , SG6043
%test 2: NREL S823 , NACA 4415 , SG6043
%test 3: NACA 4415 , SG6043    , SG6043
%test 4: NACA 4415 , NACA 4415 , SG6043

test1 = data_test(3:(end-4),[1,2,3]);


%aerofoil data
foildat = readmatrix("Aerofoildat.xlsx");
foildat(isnan(foildat)) = 0;

aerofoil_data.S823.cl = foildat(1:106,10);
aerofoil_data.SG6043.cl = foildat(1:106,18);

aerofoil_data.S823.alpha = foildat(1:106,9);
aerofoil_data.SG6043.alpha= foildat(1:106,17);

aerofoil_data.S823.cd = foildat(1:106,11);
aerofoil_data.SG6043.cd = foildat(1:106,19);

N = 200;
TSR = 4;
R = linspace(0.05,0.25,N)';
T1.R = R;

%test 1: NREL S823 , SG6043    , SG6043
%test 2: NREL S823 , NACA 4415 , SG6043
%test 3: NACA 4415 , SG6043    , SG6043
%test 4: NACA 4415 , NACA 4415 , SG6043

T1.radialstat = [0.48,0.6,1];
T1.aoa_d = deg2rad([8*ones(0.48*length(R),1);5.5*ones(0.52*length(R),1)]);
T1.cl_max = [1.1454*ones(0.48*length(R),1);1.2878*ones(0.52*length(R),1)];
T1.chord = interp1(test1(:,1),test1(:,2),R,"linear","extrap");
T1.twist = deg2rad(interp1(test1(:,1),test1(:,3),R,"linear","extrap"));
T1.aoa_dist = deg2rad([aerofoil_data.S823.alpha,aerofoil_data.SG6043.alpha,aerofoil_data.SG6043.alpha]);
T1.Cl_dist = [aerofoil_data.S823.cl,aerofoil_data.SG6043.cl,aerofoil_data.SG6043.cl];
T1.Cd_dist = [aerofoil_data.S823.cd,aerofoil_data.SG6043.cd,aerofoil_data.SG6043.cd];
% AIF1 = axial_induction_factor(TSR,T1);
% T1.a = AIF1(:,1);
% T1.cl = AIF1(:,2);
% T1.cd = AIF1(:,3);
% T1.a_prime = AIF1(:,4);
% T1.Phi = AIF1(:,5);


TSR_tested = (1:1:9)';

for i = 1:length(TSR_tested)
    AIF = axial_induction_factor(TSR,T1);
    T1.a = AIF(:,1);
    T1.cl = AIF(:,2);
    T1.cd = AIF(:,3);
    T1.a_prime = AIF(:,4);
    T1.Phi = AIF(:,5);
    
    ar = BEM_function(TSR_tested(i),T1);
    Cp_tested(i,1) = 2*ar(1);
    P_tested(i,1) = 2*ar(2);
    F_xtested(i,1) = ar(3);
    F_thettested(i,1) = ar(4);
    
end

disp(F_xtested)
disp(F_thettested)

figure;
hold on 
yline(0.593);
plot(TSR_tested,Cp_tested(:,1),Color=[0.8500 0.3250 0.0980],LineStyle="-");
title("C_P vs TSR, Test 1")
xlabel('TSR')
ylabel('C_P')
legend('betz','Test 1');
axis([0 10 0 0.65]);
hold off

figure
hold on 
yline(71.316607605748914)
plot(TSR_tested,P_tested(:,1),Color=[0.8500 0.3250 0.0980],LineStyle="-");
title("Power vs TSR, Test 1")
xlabel('TSR')
ylabel('Power (W)')
legend('betz','Test 1');
axis([0,10,0,75])
hold off

figure;
hold on 
plot(R,T1.chord)
legend('Test 1');
xlabel('Radial position (m)')
ylabel('Chord length (m)')
hold off

figure;
hold on 
plot(R,rad2deg(T1.twist))
legend('Test 1');
xlabel('Radial position (m)')
ylabel('Twist (degrees)')
hold off