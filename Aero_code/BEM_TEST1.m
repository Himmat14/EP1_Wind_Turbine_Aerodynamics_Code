clear
clc
format long

data_test = readmatrix("testdat.xlsx");

%test 1: NREL S823 , SG6043 , SG6043
%test 2: NREL S823 , NACA 4415 , SG6043
%test 3: NACA 4415 , SG6043 , SG6043
%test 4: NACA 4415 , NACA 4415 , SG6043

test1 = data_test(1:end-1,[1,2,3]);
test2 = data_test(1:end,[5,6,7]);
test3 = data_test(1:end-1,[9,10,11]);
test4 = data_test(1:end-1,[13,14,15]);

%aerofoil data
foildat = readmatrix("Aerofoildat.xlsx");
foildat(isnan(foildat)) = 0;

aerofoil_data.S823.cl = foildat(1:106,10);
aerofoil_data.SG6043.cl = foildat(1:106,18); 
aerofoil_data.NACA4415.cl = foildat(1:106,2);

aerofoil_data.S823.alpha = foildat(1:106,9);
aerofoil_data.SG6043.alpha= foildat(1:106,17); 
aerofoil_data.NACA4415.alpha = foildat(1:106,1);

aerofoil_data.S823.cd = foildat(1:106,11);
aerofoil_data.SG6043.cd = foildat(1:106,19); 
aerofoil_data.NACA4415.cd = foildat(1:106,3);

N = 200;
TSR = 3.8;
R = linspace(0.025,0.25,N)';
T1.R = R;
T2.R = R;
T3.R = R;
T4.R = R;

%test 1: NREL S823 , SG6043 , SG6043
%test 2: NREL S823 , NACA 4415 , SG6043
%test 3: NACA 4415 , SG6043 , SG6043
%test 4: NACA 4415 , NACA 4415 , SG6043

% storing + interpolating test values
T1.radialstat = [0.2,0.6,1]; 
T1.aoa_d = deg2rad([8*ones(0.2*length(R),1);5.5*ones(0.8*length(R),1)]);
T1.cl_max = [1.1454*ones(0.2*length(R),1);1.2878*ones(0.8*length(R),1)];
T1.chord = interp1(test1(:,1),test1(:,2),R,"linear","extrap");
T1.twist = deg2rad(interp1(test1(:,1),test1(:,3),R,"linear","extrap"));
T1.aoa_dist = deg2rad([aerofoil_data.S823.alpha,aerofoil_data.SG6043.alpha,aerofoil_data.SG6043.alpha]);
T1.Cl_dist = [aerofoil_data.S823.cl,aerofoil_data.SG6043.cl,aerofoil_data.SG6043.cl];
T1.Cd_dist = [aerofoil_data.S823.cd,aerofoil_data.SG6043.cd,aerofoil_data.SG6043.cd];
AIF1 = axial_induction_factor(TSR,T1);
T1.a = AIF1(:,1);
T1.cl = AIF1(:,2);
T1.cd = AIF1(:,3);
T1.a_prime = AIF1(:,4);
T1.Phi = AIF1(:,5);

T2.radialstat = [0.2,0.6,1];
T2.aoa_d = deg2rad([8*ones(0.2*length(R),1);6.5*ones(0.4*length(R),1);5.5*ones(0.4*length(R),1)]);
T2.cl_max = [1.1454*ones(0.2*length(R),1);1.1343*ones(0.4*length(R),1);1.2878*ones(0.4*length(R),1)];
T2.chord = interp1(test2(:,1),test2(:,2),R,"linear","extrap");
T2.twist = deg2rad(interp1(test2(:,1),test2(:,3),R,"linear","extrap"));
T2.aoa_dist = deg2rad([aerofoil_data.S823.alpha,aerofoil_data.NACA4415.alpha,aerofoil_data.SG6043.alpha]);
T2.Cl_dist = [aerofoil_data.S823.cl,aerofoil_data.NACA4415.cl,aerofoil_data.SG6043.cl];
T2.Cd_dist = [aerofoil_data.S823.cd,aerofoil_data.NACA4415.cd,aerofoil_data.SG6043.cd];
AIF2 = axial_induction_factor(TSR,T2);
T2.a = AIF2(:,1);
T2.cl = AIF2(:,2);
T2.cd = AIF2(:,3);
T2.a_prime = AIF2(:,4);
T2.Phi = AIF2(:,5);

%test 1: NREL S823 , SG6043 , SG6043
%test 2: NREL S823 , NACA 4415 , SG6043
%test 3: NACA 4415 , SG6043 , SG6043
%test 4: NACA 4415 , NACA 4415 , SG6043

T3.radialstat = [0.2,0.6,1];
T3.aoa_d = deg2rad([6.5*ones(0.2*length(R),1);5.5*ones(0.8*length(R),1)]);
T3.cl_max = [1.1343*ones(0.2*length(R),1);1.2878*ones(0.8*length(R),1)];
T3.chord = interp1(test3(:,1),test3(:,2),R,"linear","extrap");
T3.twist = deg2rad(interp1(test3(:,1),test3(:,3),R,"linear","extrap"));
T3.aoa_dist = deg2rad([aerofoil_data.NACA4415.alpha,aerofoil_data.SG6043.alpha,aerofoil_data.SG6043.alpha]);
T3.Cl_dist = [aerofoil_data.NACA4415.cl,aerofoil_data.SG6043.cl,aerofoil_data.SG6043.cl];
T3.Cd_dist = [aerofoil_data.NACA4415.cd,aerofoil_data.SG6043.cd,aerofoil_data.SG6043.cd];
AIF3 = axial_induction_factor(TSR,T3);
T3.a = AIF3(:,1);
T3.a(1) = 0.1;
T3.cl = AIF3(:,2);
T3.cd = AIF3(:,3);
T3.a_prime = AIF3(:,4);
T3.Phi = AIF3(:,5);

T4.radialstat = [0.2,0.6,1];
T4.aoa_d = deg2rad([6.5*ones(0.2*length(R),1);6.5*ones(0.4*length(R),1);5.5*ones(0.4*length(R),1)]);
T4.cl_max = [1.1343*ones(0.2*length(R),1);1.1343*ones(0.4*length(R),1);1.2878*ones(0.4*length(R),1)];
T4.chord = interp1(test4(:,1),test4(:,2),R,"linear","extrap");
T4.twist = deg2rad(interp1(test4(:,1),test4(:,3),R,"linear","extrap"));
T4.aoa_dist = deg2rad([aerofoil_data.NACA4415.alpha,aerofoil_data.NACA4415.alpha,aerofoil_data.SG6043.alpha]);
T4.Cl_dist = [aerofoil_data.NACA4415.cl,aerofoil_data.NACA4415.cl,aerofoil_data.SG6043.cl];
T4.Cd_dist = [aerofoil_data.NACA4415.cd,aerofoil_data.NACA4415.cd,aerofoil_data.SG6043.cd];
AIF4 = axial_induction_factor(TSR,T4);
T4.a = AIF4(:,1);
T4.cl = AIF4(:,2);
T4.cd = AIF4(:,3);
T4.a_prime = AIF4(:,4);
T4.Phi = AIF4(:,5);


%Actual data 

% arr1 = BEM_function(TSR,T1);
% disp(arr1)
% 
% arr2 = BEM_function(TSR,T2);
% disp(arr2)
% 
% arr3 = BEM_function(TSR,T3);
% disp(arr3)
% 
% arr4 = BEM_function(TSR,T4);
% disp(arr4)

tests = [T1,T2,T3,T4];
TSR_tested = (0:0.01:10)';
for j = 1
    for i = 1:length(TSR_tested)
        AIF = axial_induction_factor(TSR,tests(j));
        tests(j).a = AIF(:,1);
        tests(j).cl = AIF(:,2);
        tests(j).cd = AIF(:,3);
        tests(j).a_prime = AIF(:,4);
        tests(j).Phi = AIF(:,5);
        % tests(j).chord = AIF3(:,6);
        % tests(j).twist = AIF3(:,7); 

        ar = BEM_function(TSR_tested(i),tests(j));

        Cp_tested(i,j) = ar(1);
        P_tested(i,j) = ar(2);
        F_xtested(i,j) = ar(3);
    end
end


figure;
hold on 
yline(0.593);
plot(TSR_tested,Cp_tested(:,1),Color=[0.8500 0.3250 0.0980],LineStyle="-");
% plot(TSR_tested,Cp_tested(:,2),color=[0 0.4470 0.7410],LineStyle="--");
% plot(TSR_tested,Cp_tested(:,3),color=[0.4660 0.6740 0.1880],LineStyle=":");
% plot(TSR_tested,Cp_tested(:,4),Color='k',LineStyle='-.');
title("C_P vs TSR for different blade designs")
xlabel('TSR')
ylabel('C_P')
legend('betz','Test 1','Test 2','Test 3','Test 4');
axis([0 10 0 0.65]);
hold off


figure
hold on 
yline(71.316607605748914)
plot(TSR_tested,P_tested(:,1),Color=[0.8500 0.3250 0.0980],LineStyle="-");
% plot(TSR_tested,P_tested(:,2),color=[0 0.4470 0.7410],LineStyle="--");
% plot(TSR_tested,P_tested(:,3),color=[0.4660 0.6740 0.1880],LineStyle=":");
% plot(TSR_tested,P_tested(:,4),Color='k',LineStyle='-.');
title("Power vs TSR for different blade designs")
xlabel('TSR')
ylabel('Power (W)')
legend('betz','Test 1','Test 2','Test 3','Test 4');
axis([0,10,0,75])
hold off

%plotting

figure;
hold on 
plot(R,T1.chord)
% plot(R,T2.chord)
% plot(R,T3.chord)
% plot(R,T4.chord)
legend('Test 1','Test 2','Test 3','Test 4');
hold off

figure;
hold on 
plot(R,T1.twist)
% plot(R,T2.twist)
% plot(R,T3.twist)
% plot(R,T4.twist)
legend('Test 1','Test 2','Test 3','Test 4');
hold off





