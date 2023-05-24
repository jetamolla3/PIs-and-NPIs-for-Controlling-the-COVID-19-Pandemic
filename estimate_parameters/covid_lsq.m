function Y=covid_lsq(v,t,flags,Y0)
run_params;
if nargin>3
    y0 = Y0;
end

if ~isfield(flags,'cases')
    flags.cases=0;
end

switch flags.cases
    case 0
        params.Kc = v(1);
        params.Mc = v(2);
        params.rho0 = v(3);
        %params.rhoV0 = v(4); 
        params.p = v(4);
        %%%%%%%%%%%%%%
        params.K0 = 4*params.Kc;
        params.M0 = 2*params.Mc;
        params.rhoI = 4*params.rho0;
        params.rhoV0 = 0.5*params.rho0;%0.55
        params.rhoVI = 4*params.rhoV0;
        
end

%params.beta = params.f *(params.R0*params.phi*params.gamma/(params.alpha*(params.gamma+params.phi)+(1-params.alpha)*params.q*params.phi));
params.beta = params.f*(2*params.R0*params.phi*params.gammas*params.gammaa)/(params.gammas*params.gammaa+2*params.q*params.phi*params.gammaa+params.phi*params.gammas*(1-params.q));
odeopts = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-9);

[T,y] = ode23t(@ODEf,t,y0,odeopts,params,flags);

PM= y(:,10);
IsM = y(:,14);
IaM = y(:,18);
VPM = y(:,38);%%PM
VIsM = y(:,42);%IsM
VIaM = y(:,46);%IaM
WPM = y(:,64);%PM
WIsM = y(:,68);%IsM
WIaM = y(:,72);%IaM
M = y(:,27);


AM= PM+ IsM + IaM+ VPM+ VIsM+ VIaM+ WPM+ WIsM+ WIaM;
Vac = y(:,29)+ y(:,30)+ y(:,31)+ y(:,32)+ y(:,33)+ y(:,34)+ y(:,35)+ y(:,36)+ y(:,37)+ y(:,38)+ y(:,39)+ y(:,40)+ y(:,41)+ y(:,42)+ y(:,43)+ y(:,44)+y(:,45)+ y(:,46)+ y(:,47)+ y(:,48)+ y(:,49)+ y(:,50)+ y(:,51)+ y(:,52)+ y(:,53)+ y(:,54)+ y(:,55)+ y(:,56)+ y(:,57)+ y(:,58)+ y(:,59)+ y(:,60)+y(:,61)+ y(:,62)+ y(:,63)+ y(:,64)+ y(:,65)+ y(:,66)+ y(:,67)+ y(:,68)+ y(:,69)+ y(:,70)
%Vac = y(:,81);
%Vac = sum(y(:,29:80),2);

Y = [AM M Vac];
%Y = [AM Vac];
% Y = M;
end