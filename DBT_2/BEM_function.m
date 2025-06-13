function [parameters] = BEM_function(TSR,test_data)
    format long
    %cl_des = list of cl of each airofoil section
    %cd_des = list of cd of each aerofoil section 
    %twist = twist distribution of blade
    %chord = chord dist of blade
    %a_dist = distribution of axial induction factor across blade

    %imports all data from structure 
    CL = test_data.cl;
    CD = test_data.cd;
    Twist = test_data.twist;
    Chord = test_data.chord;
    a = test_data.a;
    R = test_data.R;
    a_prime = test_data.a_prime;
    phi = test_data.Phi;
 
    %calculates important parameters
    dR=R(2)-R(1);
    V = 10;
    rho = 1.225;
    omega = (3000/60);

    
    lambda_r = TSR.*R./R(end);
    dlambda_r = lambda_r(2)-lambda_r(1);
    P_flow = (0.5*pi*rho*(R(end).^2)*(V.^3));


    % a_prime = (V^2).*a.*(1-a)./(lambda_r.^2);
    % 
    % W = sqrt(((1-a).^2).*V^2+((R.*omega).^2).*(1+a_prime.^2));
    % 
    % mu = R./R(end);
    % 
    % dmu = mu(2)-mu(1);
    % 
    % I = sum((mu.^2).*(8.*a_prime.*(1-a).*mu - (W.*Chord.*CD.*(1+a_prime))./(V.*pi.*R)).*dmu);
    % 
    % Q = 0.5*rho*(V^2)*pi*(R(end)^3)*TSR*I;
    % 
    % P = omega.*Q/100;
    %
    % Cp = P/P_flow;
    
    F = ((2.*a.*R./pi).^2).*(cos(exp((R(end)-R)./(R(end).*sin(phi)))).*(cos(exp((R-R(1))./(R(1).*sin(phi))))));

    A = (1-(CL./CD).*cot(phi)).*((sin(phi).*lambda_r).^2);

    I2 = sum(F.*A.*(cos(phi)-lambda_r.*sin(phi)).*(sin(phi)+lambda_r.*cos(phi)).*dlambda_r);

    Cp2 = (8/(TSR^2))*I2*50;

    P2 = Cp2*P_flow;

    W = sqrt(((1-a).^2).*V^2+((R.*omega).^2).*(1+a_prime.^2));

    dF_x = rho.*(W.^2).*(CL.*sin(phi)+CD.*cos(phi)).*Chord.*dR;
    dF_theta = rho.*(W.^2).*(CL.*cos(phi)-CD.*sin(phi)).*Chord.*dR;
    F_x = sum(dF_x);
    F_theta = sum(dF_theta);

    % W = V.*(1-a)./cos(Twist);
    % 
    % dF_x = rho.*(W.^2).*(CL.*sin(phi)+CD.*cos(phi)).*Chord.*dR;
    % dF_theta = rho.*(W.^2).*(CL.*cos(phi)-CD.*sin(phi)).*Chord.*dR;
    % 
    % dT = dF_theta.*R;
    % 
    % dP = omega.*dT;
    % 
    % F_x = sum(dF_x);
    % 
    % P = sum(dP.*dR);
    % 
    % Cp = P/(0.5.*rho.*pi.*(R(end)^2).*(V.^3));

    % parameters = [Cp,P,F_x];

    % W = sqrt((omega.*R).^2+(V.*(1-a)).^2);
    % 
    % L = sum(CL.*0.5.*rho.*(W.^2).*Chord.*dR);
    % D = sum(CD.*0.5.*rho.*(W.^2).*Chord.*dR);
    % 
    % 
    % dF_x = sigma.*pi.*rho.*((W.*(1-a)./cos(Twist)).^2).*(CL.*sin(Twist)+CD.*cos(Twist)).*Chord.*R.*dR;
    % 
    % dT = sigma.*pi.*rho.*((W.*(1-a)./cos(Twist)).^2).*(CL.*cos(Twist)-CD.*sin(Twist)).*Chord.*(R.^2).*dR;

    % dF_x = 0.5*rho*(V^2).*(4.*a.*(1-a)).*2.*pi.*R.*dR;
    % dT = 4.*a_prime.*(1-a).*rho.*V.*omega.*R.*pi.*dR;

    % lambda_r = TSR.*R;
    % dlambda = lambda_r(2)-lambda_r(1);
    % 
    % Q = (2/pi).*acos(exp(-(1-R./R(end))/(R/R(end).*cos(Twist))));
    % 
    % dCp2 = Q.*(8./TSR.^2).*(lambda_r.^3).*a_prime.*(1-a).*(1-(CD./CL).*tan(Twist)).*dlambda; 
    % Cp2 = sum(dCp2);
    % 
    % dP = omega.*dT;
    % % 
    % % Cp = (lambda_r)
    % %tip loss correction
    % 
    
    %applying corrcetion factors
    % F_x = sum(dF_x);
    % P = sum(dP);
    % Cp = P./(0.5*pi*rho*(R(end).^2)*(V.^3));
    % 
    % P2 = Cp2.*(0.5*pi*rho)*(R(end)^2)*(V^3);
    % 
    % %F_x,L,D
    % 
    % dP3 = (2*a.*(1-a).*rho.*pi.*(R(end)^2).*V^3); 
    % P3 = sum(dP3);
    % Cp3 = P3./(0.5*pi*rho*(R(end).^2)*(V.^3));
    % 
    % f=(lambda_r.^3).*(a_prime.*(1-a).*(1-(CD./CL).*tan(Twist)))
    % 
    % funct = polyfit(lambda_r,((2/pi).*acos(exp(-(1-R./R(end))/(R/R(end).*cos(Twist))))).*(lambda_r.^3).*(a_prime.*(1-a).*(1-(CD./CL).*tan(Twist))),10);
    % 
    % Cp = (8/(TSR^2))*interg(lambda_r,funct);
    % 
    % P = Cp.*(0.5*pi*rho)*(R(end)^2)*(V^3);

    parameters = [Cp2,P2,F_x,F_theta];
    
end
