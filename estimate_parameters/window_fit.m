close all
clear variables

run_params;
flags.cost_control=0;
flags.model = 3;
flags.lsq_refine = 1; %turn on if you want to refine the parameter fitting

%  window(1) = 367;
%  window(2) = 393;
%  window(3) = 425;
%  window(4) = 457;
%  window(5) = 462;
%  window(6) = 476;
%%%%%%
 window(1) = 492;
 window(2) = 540;
 %%window(3) = 587;
 window(3) = 630;




load('DATA_pos');
load('DATA_tot');
load('DATA_T');
load('DATA_VAC');


spot = zeros(size(window));

for k=1:length(window)
    spot(k) = find(DATA_T==window(k));
end
%%%%%%%%%%%%%%%%%%%
%v0=[0.0031    0.0000    0.0038 0];%First part, all fixed
v0=[0.0096    0.0097    0.0043    0.0153];
%%%%%%%%%%%%%No Waning:
%v0=[0.0062    0.0066    0.0097    0.0146]
%%%%%%%%%%%%w=0:
%v0=[0.0131    0.0000    0.0003    0.0153]
%%%%%%%%%WI=0:
%v0=[0.0026    0.0028    0.0096    0.0150]
%%%%%%%%Wp=0:
% v0=[0.0098    0.0097    0.0041    0.0149]
%%%%%%%%WI=0;Wp=0:
% v0=[ 0.0062    0.0072    0.0097    0.0146]

LB = [0;0;0;0];
UB =[1;1;1;1];

flags.diffy0 = 0;

var_names = {'S';'S_1';'S_2';'E';'E_1';'E_2';'P';'P_1';'P_2';'P_M';'I_S';'I_{S1}';'I_{S2}';'I_{SM}';'I_A';'I_{A1}';'I_{A2}';'I_{AM}';'R_S';'R_{S1}';'R_{S2}';'R_{SM}';'R_A';'R_{A1}';'R_{A2}';'R_{AM}';'M';'C'};

options = optimoptions('lsqcurvefit','UseParallel',true);
odeopts = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-9);

flags.cases = 0;
for k=1:length(window)
    fprintf('Fitting window %i\n',k);
    if k==1
%         current_times = DATA_T(342:window(1));
%         current_pos = DATA_pos(342:window(1));
%         current_tot = DATA_tot(342:window(1));
%         current_Matrix = Matrix(342:window(1));
        %%%%%%%%%%
        current_times = DATA_T(476:window(1));
        current_pos = DATA_pos(476:window(1));
        current_tot = DATA_tot(476:window(1));
         current_Matrix = Matrix(476:window(1));
    else
        current_times = DATA_T(window(k-1):window(k));
        current_pos = DATA_pos(window(k-1):window(k));
        current_tot = DATA_tot(window(k-1):window(k));
        current_Matrix = Matrix(window(k-1):window(k));
        v0 = V{k-1};
        y0 = y{k-1}(end,:);
    end
    %if k==3
        %flags.cases = 1;

    %end
   if flags.lsq_refine
        %Update parameters to convergence
        converge_done = 0;
        converge_tol = 1E-8;
        while ~converge_done
            %[V{k},RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(@covid_lsq,v0,current_times,[current_pos/params.N_crit,current_Matrix/params.N_crit],LB,UB,options,flags,y0);
            [V{k},RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(@covid_lsq,v0,current_times,[current_pos/params.N_crit,current_tot/params.N_crit,current_Matrix/params.N_crit],LB,UB,options,flags,y0);
            if norm(V{k}-v0,2)<=converge_tol
                converge_done = 1;
            end
            v0 = V{k};

        end
    else
        %[V{k},RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(@covid_lsq,v0,current_times,[current_pos/params.N_crit,current_Matrix/params.N_crit],LB,UB,options,flags,y0);
        [V{k},RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(@covid_lsq,v0,current_times,[current_pos/params.N_crit,current_tot/params.N_crit,current_Matrix/params.N_crit],LB,UB,options,flags,y0);
   end
    
    switch flags.cases
        case 0
            params.Kc = V{k}(1);
            params.Mc = V{k}(2);
            params.rho0 = V{k}(3);
            params.p = V{k}(4);
            %params.w = V{k}(5);
           
            params.K0 = 4*params.Kc;
            params.M0 = 2*params.Mc;
            params.rhoI = 4*params.rho0;
            params.rhoV0 = 0.55*params.rho0;%0.55
            params.rhoVI = 4*params.rhoV0;
    end

    %params.beta = params.R0*params.phi*params.gamma/(params.alpha*(params.gamma+params.phi)+(1-params.alpha)*params.q*params.phi);
    params.beta = params.f*(2*params.R0*params.phi*params.gammas*params.gammaa)/(params.gammas*params.gammaa+2*params.q*params.phi*params.gammaa+params.phi*params.gammas*(1-params.q));
    [t{k},y{k}] = ode23t(@ODEf,current_times,y0,odeopts,params,flags);
    
    PM{k} = y{k}(:,10);
    IsM{k} = y{k}(:,14);
    IaM{k} = y{k}(:,18);
    VPM{k} = y{k}(:,38);%%PM
    VIsM{k} = y{k}(:,42);%IsM
    VIaM{k} = y{k}(:,46);%IaM
    WPM{k}= y{k}(:,64);%PM
    WIsM{k}= y{k}(:,68);%IsM
    WIaM{k}= y{k}(:,72);%IaM
    M{k} = y{k}(:,27);%Cumulative reported case
    S{k} = y{k}(:,1);
    S1{k} = y{k}(:,2);
    S2{k} = y{k}(:,3);

    AM{k} = PM{k} + IsM{k} + IaM{k}+VPM{k}+VIsM{k}+VIaM{k}+WPM{k}+WIsM{k}+WIaM{k};
    Vac{k} = y{k}(:,29)+ y{k}(:,30)+ y{k}(:,31)+ y{k}(:,32)+ y{k}(:,33)+ y{k}(:,34)+ y{k}(:,35)+ y{k}(:,36)+ y{k}(:,37)+ y{k}(:,38)+ y{k}(:,39)+ y{k}(:,40)+ y{k}(:,41)+ y{k}(:,42)+ y{k}(:,43)+ y{k}(:,44)+y{k}(:,45)+ y{k}(:,46)+ y{k}(:,47)+ y{k}(:,48)+ y{k}(:,49)+ y{k}(:,50)+ y{k}(:,51)+ y{k}(:,52)+ y{k}(:,53)+ y{k}(:,54)+ y{k}(:,55)+ y{k}(:,56)+ y{k}(:,57)+ y{k}(:,58)+ y{k}(:,59)+ y{k}(:,60)+y{k}(:,61)+ y{k}(:,62)+ y{k}(:,63)+ y{k}(:,64)+ y{k}(:,65)+ y{k}(:,66)+ y{k}(:,67)+ y{k}(:,68)+ y{k}(:,69)+ y{k}(:,70)+y{k}(:,71)+y{k}(:,72)+y{k}(:,73)+y{k}(:,74)+y{k}(:,75)+y{k}(:,76)+y{k}(:,77)+y{k}(:,78)+y{k}(:,79)+y{k}(:,80);
    AR{k} = y{k}(:,19)+y{k}(:,20)+y{k}(:,21)+y{k}(:,22)+y{k}(:,23)+y{k}(:,24)+y{k}(:,25)+y{k}(:,26)+...
        y{k}(:,47)+y{k}(:,48)+y{k}(:,49)+y{k}(:,50)+y{k}(:,51)+y{k}(:,52)+y{k}(:,53)+y{k}(:,54)+...
        y{k}(:,73)+y{k}(:,74)+y{k}(:,75)+y{k}(:,76)+y{k}(:,77)+y{k}(:,78)+y{k}(:,79)+y{k}(:,80);%All recovered population
    Recovered{k} = sum(y{k}(:,19:26),2)+sum(y{k}(:,47:54),2)+sum(y{k}(:,73:80));
    tests{k} = params.rho0*(sum(y{k}(:,1:10),2)+sum(y{k}(:,15:26),2))+params.rhoI*(sum(y{k}(:,11:14),2))+params.rhoV0*(sum(y{k}(:,29:38),2)+sum(y{k}(:,43:64),2)+sum(y{k}(:,69:80),2))+params.rhoVI*(sum(y{k}(:,39:42),2)+sum(y{k}(:,65:68),2));
    AR1{k} = y{k}(:,19)+y{k}(:,20)+y{k}(:,21)+y{k}(:,22)+...
        y{k}(:,47)+y{k}(:,48)+y{k}(:,49)+y{k}(:,50)+...
        y{k}(:,73)+y{k}(:,74)+y{k}(:,75)+y{k}(:,76);
    %%%
     AR2{k} = y{k}(:,23)+y{k}(:,24)+y{k}(:,25)+y{k}(:,26)+...
              y{k}(:,51)+y{k}(:,52)+y{k}(:,53)+y{k}(:,54)+...
              y{k}(:,77)+y{k}(:,78)+y{k}(:,79)+y{k}(:,80);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ascertainment:
    P{k} = y{k}(:,7);
    P1{k} = y{k}(:,8);
    P2{k} = y{k}(:,9);
    Is{k} = y{k}(:,11);
    Is1{k} = y{k}(:,12);
    Is2{k} = y{k}(:,13);
    Ia{k} = y{k}(:,15);
    Ia1{k} = y{k}(:,16);
    Ia2{k} = y{k}(:,17);
    %%%%
    VP{k} = y{k}(:,35);
    VP1{k} = y{k}(:,36);
    VP2{k} = y{k}(:,37);
    VIs{k} = y{k}(:,39);
    VIs1{k} = y{k}(:,40);
    VIs2{k} = y{k}(:,41);
    VIa{k} = y{k}(:,43);
    VIa1{k} = y{k}(:,44);
    VIa2{k} = y{k}(:,45);
    %%%%
    WP{k} = y{k}(:,61);
    WP1{k} = y{k}(:,62);
    WP2{k} = y{k}(:,63);
    WIs{k} = y{k}(:,65);
    WIs1{k} = y{k}(:,66);
    WIs2{k} = y{k}(:,67);
    WIa{k} = y{k}(:,69);
    WIa1{k} = y{k}(:,70);
    WIa2{k} = y{k}(:,71);
    In{k}=(PM{k} + IsM{k} + IaM{k}+VPM{k}+VIsM{k}+VIaM{k}+WPM{k}+WIsM{k}+WIaM{k})/(params.rho0*(Ia{k}+Ia1{k}+Ia2{k}+P{k}+P1{k}+P2{k}) + params.rhoI*(Is{k}+Is1{k}+Is2{k}) + params.rhoV0*(VIa{k}+WIa{k}+VIa1{k}+WIa1{k}+VIa2{k}+WIa2{k}+VP{k}+WP{k}+VP1{k}+WP1{k}+VP2{k}+WP2{k}) + params.rhoVI*(VIs{k}+WIs{k}+VIs1{k}+WIs1{k}+VIs2{k}+WIs2{k}));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mu function:
    %%%
    %KM{k} = min(1/log(2)*max((params.rho0*(Ia{k}+Ia1{k}+Ia2{k}+P{k}+P1{k}+P2{k}) + params.rhoI*(Is{k}+Is1{k}+Is2{k}) + params.rhoV0*(VIa{k}+WIa{k}+VIa1{k}+WIa1{k}+VIa2{k}+WIa2{k}+VP{k}+WP{k}+VP1{k}+WP1{k}+VP2{k}+WP2{k}) + params.rhoVI*(VIs{k}+WIs{k}+VIs1{k}+WIs1{k}+VIs2{k}+WIs2{k}))/M{k},0),1000);
    %mu{k} = max(KM{k}-params.Kc,0)/(max(KM{k}-params.Kc,0)+params.K0-params.Kc)*max(AM{k}-params.Mc,0)/(max(AM{k}-params.Mc,0)+params.M0-params.Mc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot(t{1},mu{1},'linewidth',2);
% hold on
% for k=2:length(window)
%     plot(t{k},mu{k},'linewidth',2);
% end
% hold off
% title('Transition fuction:mu');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(t{1},AM{1},'linewidth',2);
hold on
for k=2:length(window)
    plot(t{k},AM{k},'linewidth',2);
end
plot(DATA_T,(DATA_pos)/params.N0,'o');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(t{1},M{1},'linewidth',2);
hold on
for k=2:length(window)
    plot(t{k},M{k},'linewidth',2);
end
plot(DATA_T,(DATA_tot)/params.N0,'o');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(t{1},Vac{1},'linewidth',2);
hold on
for k=2:length(window)
    plot(t{k},Vac{k},'linewidth',2);
end
plot(DATA_T,Matrix/params.N0,'o');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(t{1},AR{1}*params.N0/params.N,'linewidth',2);
hold on
for k=2:length(window)
    plot(t{k},AR{k}*(params.N0/params.N)*100,'linewidth',2);
end
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
for k=1:length(window)
    plot(t{k}(2:end),(M{k}(2:end)-M{k}(1:end-1))./tests{k}(2:end),'linewidth',2);
end
hold off
title('Positivity Rate');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
for k=1:length(window)
    plot(t{k}(2:end),(M{k}(2:end)-M{k}(1:end-1)),'linewidth',2);
end
hold off
title('Daily Incidence');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%y{3}(:,2);y{3}(20,2)