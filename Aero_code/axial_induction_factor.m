%itteration function
function [parameters] = axial_induction_factor(TSR,test_data)

    Twist = test_data.twist;
    Chord = test_data.chord;
    R = test_data.R;
    AOA = test_data.aoa_dist;
    Cl = test_data.Cl_dist;
    Cd = test_data.Cd_dist;
    station = test_data.radialstat;
    aoa_dist = test_data.aoa_d;
    cl_des = test_data.cl_max;

    a = linspace(1/3,1/3,length(R))';

    lambda_r = TSR.*R./(R(end));
    sigma_r = Chord./(pi.*R);
    V = 10;

    theta = zeros(length(R),1);
    alpha = zeros(length(R),1);
    CL = zeros(length(R),1);
    CD = zeros(length(R),1);
    
    

    
    for x = 1:10000
        for i = 1:(station(1)*length(R))
            theta(i) = atan((1-a(i))/lambda_r(i));
            alpha(i) = theta(i)-Twist(i);
            cl_local = interp1(AOA(:,1),Cl(:,1),alpha(i),"linear","extrap");
            CL(i) = cl_local;
            CD(i) = interp1(AOA(:,1),Cd(:,1),alpha(i),"linear","extrap");
            a(i) = (sigma_r(i)*CL(i)/(4*sin(theta(i)*tan(theta(i)))))*(1-a(i));
            theta2(i) = atan((1-a(i))/lambda_r(i));
        end
        
        for i = (station(1)*length(R)):(station(2)*length(R))
            theta(i) = atan((1-a(i))/lambda_r(i));
            alpha(i) = theta(i)-Twist(i);
            cl_local = interp1(AOA(:,2),Cl(:,2),alpha(i),"linear","extrap");
            CL(i) = cl_local;
            CD(i) = interp1(AOA(:,2),Cd(:,2),alpha(i),"linear","extrap");
            a(i) = (sigma_r(i)*CL(i)/(4*sin(theta(i)*tan(theta(i)))))*(1-a(i));
            theta2(i) = atan((1-a(i))/lambda_r(i));
        end

        for i = (station(2)*length(R)):(station(3)*length(R))
            theta(i) = atan((1-a(i))/lambda_r(i));
            alpha(i) = theta(i)-Twist(i);
            cl_local = interp1(AOA(:,3),Cl(:,3),alpha(i),"linear","extrap");
            CL(i) = cl_local;
            CD(i) = interp1(AOA(:,3),Cd(:,3),alpha(i),"linear","extrap");
            a(i) = (sigma_r(i)*CL(i)/(4*sin(theta(i)*tan(theta(i)))))*(1-a(i));
            theta2(i) = atan((1-a(i))/lambda_r(i));
        end
    end


    
    phi = theta;
    % twi = theta - aoa_dist;
     
    % chord = (pi.*R.*8)/(9.*cl_des.*TSR.*sqrt((TSR^2).*(R./R(end)).^2)+4/9);
    a_prime = (V^2).*a.*(1-a)./(lambda_r.^2);
     
    parameters = [a,CL,CD,a_prime,phi];

end

