function dy=ODEf(t,y,p,flags)

switch flags.model
        
    case 1

        
    case 2

        
    case 3
        dy = zeros(80,1);
        
        S = y(1);
        S1 = y(2);
        S2 = y(3);
        E = y(4);
        E1 = y(5);
        E2 = y(6);
        P = y(7);
        P1 = y(8);
        P2 = y(9);
        PM = y(10);
        Is = y(11);
        Is1 = y(12);
        Is2 = y(13);
        IsM = y(14);
        Ia = y(15);
        Ia1 = y(16);
        Ia2 = y(17);
        IaM = y(18);
        Rs = y(19);
        Rs1 = y(20);
        Rs2 = y(21);
        RsM = y(22);
        Ra = y(23);
        Ra1 = y(24);
        Ra2 = y(25);
        RaM = y(26);
        M = y(27);
        C = y(28);
        VS = y(29);
        VS1 = y(30);
        VS2 = y(31);
        VE = y(32);
        VE1 = y(33);
        VE2 = y(34);
        VP = y(35);
        VP1 = y(36);
        VP2 = y(37);
        VPM = y(38);
        VIs = y(39);
        VIs1 = y(40);
        VIs2 = y(41);
        VIsM = y(42);
        VIa = y(43);
        VIa1 = y(44);
        VIa2 = y(45);
        VIaM = y(46);
        VRs = y(47);
        VRs1 = y(48);
        VRs2 = y(49);
        VRsM = y(50);
        VRa = y(51);
        VRa1 = y(52);
        VRa2 = y(53);
        VRaM = y(54);
        WS = y(55);
        WS1 = y(56);
        WS2 = y(57);
        WE = y(58);
        WE1 = y(59);
        WE2 = y(60);
        WP = y(61);
        WP1 = y(62);
        WP2 = y(63);
        WPM = y(64);
        WIs = y(65);
        WIs1 = y(66);
        WIs2 = y(67);
        WIsM = y(68);
        WIa = y(69);
        WIa1 = y(70);
        WIa2 = y(71);
        WIaM = y(72);
        WRs = y(73);
        WRs1 = y(74);
        WRs2 = y(75);
        WRsM = y(76);
        WRa = y(77);
        WRa1 = y(78);
        WRa2 = y(79);
        WRaM = y(80);
        
        
        Infectious = p.alpha*(P+WP)+p.alpha*(Ia+WIa)+(Is+WIs)...
            +p.alpha*p.delta*(P1+WP1)+p.alpha*p.delta*(Ia1+WIa1)+p.delta*(Is1+WIs1)...
            +p.zeta*(p.alpha*VP+p.alpha*VIa+VIs+p.alpha*p.delta*VP1...
            +p.alpha*p.delta*VIa1+p.delta*VIs1);
        
        F_S = p.N0/p.N*p.beta*S*Infectious;
        F_S1 = p.N0/p.N*p.beta*p.delta*S1*Infectious;
        F_VS = p.N0/p.N*p.beta*p.epsilon*VS*Infectious;
        F_VS1 = p.N0/p.N*p.beta*p.epsilon*p.delta*VS1*Infectious;
        F_WS = p.N0/p.N*p.beta*WS*Infectious;
        F_WS1 = p.N0/p.N*p.beta*p.delta*WS1*Infectious;
        
            KM = min(1/log(2)*max((p.rho0*(Ia+Ia1+Ia2+P+P1+P2) + p.rhoI*(Is+Is1+Is2) + p.rhoV0*(VIa+WIa+VIa1+WIa1+VIa2+WIa2+VP+WP+VP1+WP1+VP2+WP2) + p.rhoVI*(VIs+WIs+VIs1+WIs1+VIs2+WIs2))/M,0),1000);
            AM = PM+IsM+IaM+VPM+VIsM+WIaM+WPM+WIsM+WIaM;
%             AM = p.rho0*(P+P1+P2)+p.rhoI*(Is+Is1+Is2)+p.rho0*(Ia+Ia1+Ia2);

            mu = max(KM-p.Kc,0)/(max(KM-p.Kc,0)+p.K0-p.Kc)*max(AM-p.Mc,0)/(max(AM-p.Mc,0)+p.M0-p.Mc);
            
            
            if flags.cost_control
                nu = max(C-p.Cc,0)/(max(C-p.Cc,0)+p.C0-p.Cc)*max(p.eta*p.Mc-AM,0);%/p.eta/p.Mc;
            else
                nu = max(C-p.Cc,0)/(max(C-p.Cc,0)+p.C0-p.Cc);
            end
            
        

        
        dy(1) = -F_S - p.mumax*mu*S + p.numax/2*nu*S1 + (1-p.q2)*p.numax*nu*S2 - p.p*S+p.w*WS+p.wI*(Rs+Ra);
        dy(2) = -F_S1 - p.mumax/2*mu*S1 + p.q1*p.mumax*mu*S - p.numax/2*nu*S1 + p.q2*p.numax*nu*S2 - p.p*S1+p.w*WS1+p.wI*(Rs1+Ra1);
        dy(3) = (1-p.q1)*p.mumax*mu*S + p.mumax/2*mu*S1 - p.numax*nu*S2 - p.p*S2+p.w*WS2+p.wI*(Rs2+Ra2);
        dy(4) = F_S - p.mumax*mu*E + p.numax/2*nu*E1 + (1-p.q2)*p.numax*nu*E2 -p.sigma*E - p.p*E+p.w*WE;
        dy(5) = F_S1 - p.mumax/2*mu*E1 + p.q1*p.mumax*mu*E - p.numax/2*nu*E1 + p.q2*p.numax*nu*E2 - p.sigma*E1 - p.p*E1+p.w*WE1;
        dy(6) = (1-p.q1)*p.mumax*mu*E + p.mumax/2*mu*E1 - p.numax*nu*E2 - p.sigma*E2 - p.p*E2+p.w*WE2;
        dy(7) = p.sigma*E -p.mumax*mu*P + p.numax/2*nu*P1 + (1-p.q2)*p.numax*nu*P2 - p.phi*P - p.rho0*P - p.p*P+p.w*WP;
        dy(8) = p.sigma*E1 - p.mumax/2*mu*P1 + p.q1*p.mumax*mu*P - p.numax/2*nu*P1 + p.q2*p.numax*nu*P2 - p.phi*P1 - p.rho0*P1 - p.p*P1+p.w*WP1;
        dy(9) = p.sigma*E2 + (1-p.q1)*p.mumax*mu*P + p.mumax/2*mu*P1 - p.numax*nu*P2 - p.phi*P2 - p.rho0*P2 - p.p*P2+p.w*WP2;
        dy(10) = p.rho0*(P+P1+P2) - p.phi*PM+p.w*WPM;
        dy(11) = p.q*p.phi*P - p.muI*Is - p.gammas*Is - p.rhoI*Is+p.w*WIs;
        dy(12) = p.q*p.phi*P1 + p.qI*p.muI*Is - p.gammas*Is1 - p.rhoI*Is1+p.w*WIs1;
        dy(13) = p.q*p.phi*P2 + (1-p.qI)*p.muI*Is - p.gammas*Is2 - p.rhoI*Is2+p.w*WIs2;
        dy(14) = p.rhoI*(Is+Is1+Is2) + p.q*p.phi*PM - p.gammas*IsM+p.w*WIsM;
        dy(15) = (1-p.q)*p.phi*P - p.mumax*mu*Ia + p.numax/2*nu*Ia1 + (1-p.q2)*p.numax*nu*Ia2 - p.gammaa*Ia - p.rho0*Ia - p.p*Ia+p.w*WIa;
        dy(16) = (1-p.q)*p.phi*P1 - p.mumax/2*mu*Ia1 + p.q1*p.mumax*mu*Ia - p.numax/2*nu*Ia1 + p.q2*p.numax*nu*Ia2 - p.gammaa*Ia1 - p.rho0*Ia1 - p.p*Ia1+p.w*WIa1;
        dy(17) = (1-p.q)*p.phi*P2 + (1-p.q1)*p.mumax*mu*Ia + p.mumax/2*mu*Ia1 - p.numax*nu*Ia2 - p.gammaa*Ia2 - p.rho0*Ia2 - p.p*Ia2+p.w*WIa2;
        dy(18) = p.rho0*(Ia+Ia1+Ia2) +(1-p.q)*p.phi*PM - p.gammaa*IaM+p.w*WIaM;
        dy(19) = p.gammas*Is - p.p*Rs+p.w*WRs-p.wI*Rs- p.mumax*mu*Rs + p.numax/2*nu*Rs1 + (1-p.q2)*p.numax*nu*Rs2;
        dy(20) = p.gammas*Is1 - p.p*Rs1+p.w*WRs1-p.wI*Rs1- p.mumax/2*mu*Rs1 + p.q1*p.mumax*mu*Rs - p.numax/2*nu*Rs1 + p.q2*p.numax*nu*Rs2;
        dy(21) = p.gammas*Is2 - p.p*Rs2+p.w*WRs2-p.wI*Rs2 +(1-p.q1)*p.mumax*mu*Rs + p.mumax/2*mu*Rs1 - p.numax*nu*Rs2;
        dy(22) = p.gammas*IsM - p.p*RsM+p.w*WRsM-p.wI*RsM;
        dy(23) = p.gammaa*Ia - p.p*Ra+p.w*WRa-p.wI*Ra- p.mumax*mu*Ra + p.numax/2*nu*Ra1 + (1-p.q2)*p.numax*nu*Ra2;
        dy(24) = p.gammaa*Ia1 - p.p*Ra1+p.w*WRa1-p.wI*Ra1 - p.mumax/2*mu*Ra1 + p.q1*p.mumax*mu*Ra - p.numax/2*nu*Ra1 + p.q2*p.numax*nu*Ra2;
        dy(25) = p.gammaa*Ia2 - p.p*Ra2+p.w*WRa2-p.wI*Ra2+(1-p.q1)*p.mumax*mu*Ra + p.mumax/2*mu*Ra1 - p.numax*nu*Ra2;
        dy(26) = p.gammaa*IaM - p.p*RaM+p.w*WRaM-p.wI*RaM; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dy(27) = p.rho0*(Ia+Ia1+Ia2+P+P1+P2) + p.rhoI*(Is+Is1+Is2) + p.rhoV0*(VIa+WIa+VIa1+WIa1+VIa2+WIa2+VP+WP+VP1+WP1+VP2+WP2) + p.rhoVI*(VIs+WIs+VIs1+WIs1+VIs2+WIs2);
        dy(28) = p.n*((S2+E2+VS2+VE2+WS2+WE2)+(1-p.delta)*(S1+E1+VS1+VE1+WE1+WE2)) - p.muC*C;
        
        %Vaccination compartments start here
        dy(29) = -F_VS - p.mumaxV*mu*VS + p.numaxV/2*nu*VS1 + (1-p.q2)*p.numaxV*nu*VS2 + p.qVS*p.p*S-p.wp*VS;
        dy(30) = -F_VS1 - p.mumaxV/2*mu*VS1 + p.q1*p.mumaxV*mu*VS - p.numaxV/2*nu*VS1 + p.q2*p.numaxV*nu*VS2 + p.qVS*p.p*S1-p.wp*VS1;
        dy(31) = (1-p.q1)*p.mumaxV*mu*VS + p.mumaxV/2*mu*VS1 - p.numaxV*nu*VS2 + p.qVS*p.p*S2-p.wp*VS2;
        dy(32) = F_VS - p.mumaxV*mu*VE + p.numaxV/2*nu*VE1 + (1-p.q2)*p.numaxV*nu*VE2 -p.sigma*VE + p.qVE*p.p*E;
        dy(33) = F_VS1 - p.mumaxV/2*mu*VE1 + p.q1*p.mumaxV*mu*VE - p.numaxV/2*nu*VE1 + p.q2*p.numaxV*nu*VE2 - p.sigma*VE1 + p.qVE*p.p*E1;
        dy(34) = (1-p.q1)*p.mumaxV*mu*VE + p.mumaxV/2*mu*VE1 - p.numaxV*nu*VE2 - p.sigma*VE2 + p.qVE*p.p*E2;
        dy(35) = p.sigma*VE -p.mumaxV*mu*VP + p.numaxV/2*nu*VP1 + (1-p.q2)*p.numaxV*nu*VP2 - p.phi*VP - p.rhoV0*VP + p.qVP*p.p*P;
        dy(36) = p.sigma*VE1 - p.mumaxV/2*mu*VP1 + p.q1*p.mumaxV*mu*VP - p.numaxV/2*nu*VP1 + p.q2*p.numaxV*nu*VP2 - p.phi*VP1 - p.rhoV0*VP1 + p.qVP*p.p*P1;
        dy(37) = p.sigma*VE2 + (1-p.q1)*p.mumaxV*mu*VP + p.mumaxV/2*mu*VP1 - p.numaxV*nu*VP2 - p.phi*VP2 - p.rhoV0*VP2 + p.qVP*p.p*P2;
        dy(38) = p.rhoV0*(VP+VP1+VP2) - p.phi*VPM;
        dy(39) = p.qv*p.phi*VP - p.muI*VIs - p.gammas*VIs - p.rhoVI*VIs;
        dy(40) = p.qv*p.phi*VP1 + p.qI*p.muI*VIs - p.gammas*VIs1 - p.rhoVI*VIs1;
        dy(41) = p.qv*p.phi*VP2 + (1-p.qI)*p.muI*VIs - p.gammas*VIs2 - p.rhoVI*VIs2;
        dy(42) = p.rhoVI*(VIs+VIs1+VIs2) + p.qv*p.phi*VPM - p.gammas*VIsM;
        dy(43) = (1-p.qv)*p.phi*VP - p.mumaxV*mu*VIa + p.numaxV/2*nu*VIa1 + (1-p.q2)*p.numaxV*nu*VIa2 - p.gammaa*VIa - p.rhoV0*VIa + p.qVIa*p.p*Ia;
        dy(44) = (1-p.qv)*p.phi*VP1 - p.mumaxV/2*mu*VIa1 + p.q1*p.mumaxV*mu*VIa - p.numaxV/2*nu*VIa1 + p.q2*p.numaxV*nu*VIa2 - p.gammaa*VIa1 - p.rhoV0*VIa1 + p.qVIa*p.p*Ia1;
        dy(45) = (1-p.qv)*p.phi*VP2 + (1-p.q1)*p.mumaxV*mu*VIa + p.mumaxV/2*mu*VIa1 - p.numaxV*nu*VIa2 - p.gammaa*VIa2 - p.rhoV0*VIa2 + p.qVIa*p.p*Ia2;
        dy(46) = p.rhoV0*(VIa+VIa1+VIa2) +(1-p.qv)*p.phi*VPM - p.gammaa*VIaM;
        dy(47) = p.gammas*VIs + p.qVR*p.p*Rs-p.wI*VRs- p.mumax*mu*VRs + p.numax/2*nu*VRs1 + (1-p.q2)*p.numax*nu*VRs2;
        dy(48) = p.gammas*VIs1 + p.qVR*p.p*Rs1-p.wI*VRs1 - p.mumax/2*mu*VRs1 + p.q1*p.mumax*mu*VRs - p.numax/2*nu*VRs1 + p.q2*p.numax*nu*VRs2;
        dy(49) = p.gammas*VIs2 + p.qVR*p.p*Rs2-p.wI*VRs2 +(1-p.q1)*p.mumax*mu*VRs + p.mumax/2*mu*VRs1 - p.numax*nu*VRs2;
        dy(50) = p.gammas*VIsM + p.qVR*p.p*RsM-p.wI*VRsM;
        dy(51) = p.gammaa*VIa + p.qVR*p.p*Ra-p.wI*VRa- p.mumax*mu*VRa + p.numax/2*nu*VRa1 + (1-p.q2)*p.numax*nu*VRa2;
        dy(52) = p.gammaa*VIa1 + p.qVR*p.p*Ra1-p.wI*VRa1- p.mumax/2*mu*VRa1 + p.q1*p.mumax*mu*VRa - p.numax/2*nu*VRa1 + p.q2*p.numax*nu*VRa2;
        dy(53) = p.gammaa*VIa2 + p.qVR*p.p*Ra2-p.wI*VRa2+(1-p.q1)*p.mumax*mu*VRa + p.mumax/2*mu*VRa1 - p.numax*nu*VRa2;
        dy(54) = p.gammaa*VIaM + p.qVR*p.p*RaM-p.wI*VRaM;     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Waning:
        dy(55) = -F_WS - p.mumaxV*mu*WS + p.numaxV/2*nu*WS1 + (1-p.q2)*p.numaxV*nu*WS2 + (1-p.qVS)*p.p*S-p.w*WS+p.wI*(VRs+VRa+WRs+WRa+VRsM+VRaM+WRsM+WRaM+RsM+RaM)+p.wp*VS;
        dy(56) = -F_WS1 - p.mumaxV/2*mu*WS1 + p.q1*p.mumaxV*mu*WS - p.numaxV/2*nu*WS1 + p.q2*p.numaxV*nu*WS2 + (1-p.qVS)*p.p*S1-p.w*WS1+p.wI*(VRs1+VRa1+WRs1+WRa1)+p.wp*VS1;
        dy(57) = (1-p.q1)*p.mumaxV*mu*WS + p.mumaxV/2*mu*WS1 - p.numaxV*nu*WS2 + (1-p.qVS)*p.p*S2-p.w*WS2+p.wI*(VRs2+VRa2+WRs2+WRa2)+p.wp*VS2;
        dy(58) = F_WS - p.mumaxV*mu*WE + p.numaxV/2*nu*WE1 + (1-p.q2)*p.numaxV*nu*WE2 - p.sigma*WE + (1-p.qVE)*p.p*E-p.w*WE;
        dy(59) = F_WS1 - p.mumaxV/2*mu*WE1 + p.q1*p.mumaxV*mu*WE - p.numaxV/2*nu*WE1 + p.q2*p.numaxV*nu*WE2 - p.sigma*WE1 + (1-p.qVE)*p.p*E1-p.w*WE1;
        dy(60) = (1-p.q1)*p.mumaxV*mu*WE + p.mumaxV/2*mu*WE1 - p.numaxV*nu*WE2 - p.sigma*WE2 + (1-p.qVE)*p.p*E2-p.w*WE2;
        dy(61) = p.sigma*WE -p.mumaxV*mu*WP + p.numaxV/2*nu*WP1 + (1-p.q2)*p.numaxV*nu*WP2 - p.phi*WP - p.rhoV0*WP + (1-p.qVP)*p.p*P-p.w*WP;
        dy(62) = p.sigma*WE1 - p.mumaxV/2*mu*WP1 + p.q1*p.mumaxV*mu*WP - p.numaxV/2*nu*WP1 + p.q2*p.numaxV*nu*WP2 - p.phi*WP1 - p.rhoV0*WP1 + (1-p.qVP)*p.p*P1-p.w*WP1;
        dy(63) = p.sigma*WE2 + (1-p.q1)*p.mumaxV*mu*WP + p.mumaxV/2*mu*WP1 - p.numaxV*nu*WP2 - p.phi*WP2 - p.rhoV0*WP2 + (1-p.qVP)*p.p*P2-p.w*WP2;
        dy(64) = p.rhoV0*(WP+WP1+WP2) - p.phi*WPM-p.w*WPM;
        dy(65) = p.q*p.phi*WP - p.muI*WIs - p.gammas*WIs - p.rhoVI*WIs-p.w*WIs;
        dy(66) = p.q*p.phi*WP1 + p.qI*p.muI*WIs - p.gammas*WIs1 - p.rhoVI*WIs1-p.w*WIs1;
        dy(67) = p.q*p.phi*WP2 + (1-p.qI)*p.muI*WIs - p.gammas*WIs2 - p.rhoVI*WIs2-p.w*WIs2;
        dy(68) = p.rhoVI*(WIs+WIs1+WIs2) + p.q*p.phi*WPM - p.gammas*WIsM-p.w*WIsM;
        dy(69) = (1-p.q)*p.phi*WP - p.mumaxV*mu*WIa + p.numaxV/2*nu*WIa1 + (1-p.q2)*p.numaxV*nu*WIa2 - p.gammaa*WIa - p.rhoV0*WIa + (1-p.qVIa)*p.p*Ia-p.w*WIa;
        dy(70) = (1-p.q)*p.phi*WP1 - p.mumaxV/2*mu*WIa1 + p.q1*p.mumaxV*mu*WIa - p.numaxV/2*nu*WIa1 + p.q2*p.numaxV*nu*WIa2 - p.gammaa*WIa1 - p.rhoV0*WIa1 + (1-p.qVIa)*p.p*Ia1-p.w*WIa1;
        dy(71) = (1-p.q)*p.phi*WP2 + (1-p.q1)*p.mumaxV*mu*WIa + p.mumaxV/2*mu*WIa1 - p.numaxV*nu*WIa2 - p.gammaa*WIa2 - p.rhoV0*WIa2 + (1-p.qVIa)*p.p*Ia2-p.w*WIa2;
        dy(72) = p.rhoV0*(WIa+WIa1+WIa2) +(1-p.q)*p.phi*WPM - p.gammaa*WIaM-p.w*WIaM;
        dy(73) = p.gammas*WIs + (1-p.qVR)*p.p*Rs-p.w*WRs-p.wI*WRs - p.mumax*mu*WRs + p.numax/2*nu*WRs1 + (1-p.q2)*p.numax*nu*WRs2;
        dy(74) = p.gammas*WIs1 + (1-p.qVR)*p.p*Rs1-p.w*WRs1-p.wI*WRs1- p.mumax/2*mu*WRs1 + p.q1*p.mumax*mu*WRs - p.numax/2*nu*WRs1 + p.q2*p.numax*nu*WRs2;
        dy(75) = p.gammas*WIs2 + (1-p.qVR)*p.p*Rs2-p.w*WRs2-p.wI*WRs2+(1-p.q1)*p.mumax*mu*WRs + p.mumax/2*mu*WRs1 - p.numax*nu*WRs2;
        dy(76) = p.gammas*WIsM + (1-p.qVR)*p.p*RsM-p.w*WRsM-p.wI*WRsM;
        dy(77) = p.gammaa*WIa + (1-p.qVR)*p.p*Ra-p.w*WRa-p.wI*WRa- p.mumax*mu*WRa + p.numax/2*nu*WRa1 + (1-p.q2)*p.numax*nu*WRa2;
        dy(78) = p.gammaa*WIa1 + (1-p.qVR)*p.p*Ra1-p.w*WRa1-p.wI*WRa1- p.mumax/2*mu*WRa1 + p.q1*p.mumax*mu*WRa - p.numax/2*nu*WRa1 + p.q2*p.numax*nu*WRa2;
        dy(79) = p.gammaa*WIa2 + (1-p.qVR)*p.p*Ra2-p.w*WRa2-p.wI*WRa2+(1-p.q1)*p.mumax*mu*WRa + p.mumax/2*mu*WRa1 - p.numax*nu*WRa2;
        dy(80) = p.gammaa*WIaM + (1-p.qVR)*p.p*RaM-p.w*WRaM-p.wI*WRaM;

end

