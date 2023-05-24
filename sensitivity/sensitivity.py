import numpy as np
import pandas as pd
import logging 
from param import PARAMS 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import  sys
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

   
def  old_model(y,t,v0,v1,v2,v3):
    
    N0 = PARAMS['N0']; Kc =v0 ;Mc = v1;rho0 =v2 ; p = v3
   
    alpha = PARAMS['alpha']; sigma = PARAMS['sigma'] ; phi = PARAMS['phi'] ; gamma =PARAMS['gamma'] 
    
    q =  PARAMS['q']  ; q1 = PARAMS['q1'] ; q2 = PARAMS['q2'] ; muI = PARAMS['muI'] ; qI = PARAMS['qI'] ; qv = PARAMS['qv'] ; qw =  PARAMS['qw']
    rhoI = 4.*rho0; N_crit = PARAMS['N_crit'] ; delta = PARAMS['delta'] ; muI_v = 0.5*muI ; 
    numax = PARAMS['numax'] ; mumax = PARAMS['mumax'] ; muC = PARAMS['muC']; mumaxV = 0.5*mumax; numaxV = 2.*numax
    N = PARAMS['N'] ; n = N0/N ; T = PARAMS['time_sim']
    M0 = 2.*Mc ; K0= 4.*Kc ;  rho0_v = 0.5*rho0 ; rhoI_v = 4.*rho0_v ; rhoI_w = rhoI_v
    q_v =  PARAMS['q_v'] ; rhoV0 = rho0_v ; rhoVI = rhoI_v
    rho0_w = rho0 
    qVS = PARAMS['qVS'] ; qVE = qVS ; qVP = qVS ; qVIa = qVS ; qVR = qVS
    Cc = PARAMS['Cc'] ; C0 = 2.*Cc  ; R0 =  PARAMS['R0'] ; control = PARAMS['control'] ; 
    beta = R0*phi*gamma/(alpha*(gamma+phi)+(1-alpha)*q*phi)
    beta_s = beta*(N0/N)
    #waning rates 
   # wA = 0.005 ; wV = 0.005 ;  wI =0.005; w = 0.07 ;  wp = 0.005
    global wA ; global wp ; global w ; global epsilon ; global wV ; global wI
    
    
  
    S, S1, S2,E, E1, E2,P, P1, P2,PM,Is,Is1, Is2,IsM, Ia,Ia1,Ia2, IaM, Rs,Rs1,Rs2,RsM,Ra,Ra1,Ra2,RaM, M, C,\
    VS, VS1, VS2, VE, VE1, VE2, VP, VP1, VP2, VPM, VIs, VIs1, VIs2, VIsM, VIa,VIa1,VIa2, VIaM, VRs,VRs1,VRs2, VRsM, VRa, VRa1, VRa2, VRaM,\
    WS, WS1, WS2, WE, WE1, WE2, WP, WP1, WP2, WPM, WIs, WIs1, WIs2, WIsM, WIa,WIa1,WIa2, WIaM, WRs,WRs1,WRs2, WRsM, WRa, WRa1,WRa2, WRaM,\
    Vac, D_cases, Pos= y 
         
    def F_s(t):
        
        return (alpha*(P+WP+Ia+WIa)+(Is+WIs)+alpha*delta*(P1+ WP1 +Ia1+WIa1)+delta*(Is1+WIs1))
               
    def beta_t(t): 
        if t<=336:
            beta_f = beta_s
        elif 336<t<=643:
            beta_f = 1.5*beta_s
        else: 
            beta_f = 3.*beta_s
            
        return beta_f
            
    
    
    def beta_v(t): 
        if t<=336:
            beta_f = beta_s
        elif 336<t<=643:
            beta_f = 1.5*beta_s
        else: 
            beta_f = 3.*beta_s
            
        return beta_f
                               
    def F_sv(t):
        
        return (alpha*(VP + VIa + delta*(VP1+VIa1))+VIs +delta*VIs1)
                                                        

    def mu(t):
    
       
        KM = min(1./np.log(2)*max((rho0*(Ia+Ia1+Ia2+P+P1+P2) + rhoI*(Is+Is1+Is2 ) +rhoV0*(VIa+WIa+VIa1+WIa1+VIa2+WIa2+VP+WP+VP1+WP1+VP2+WP2) \
                                  + rhoVI*(VIs+WIs+VIs1+WIs1+VIs2+WIs2) )/M,0),1000);
            
        AM =  PM+IsM+IaM+VPM+VIsM+WIaM+WPM+WIsM+WIaM #+H +Hv

    
        
        Mu = mumax*(max(KM-Kc,0)/(max(KM-Kc,0)+K0-Kc))*(max(AM-Mc,0)/(max(AM-Mc,0)+M0-Mc))
        
    
        return Mu
    
    control = 'control'
    def nu(t):
        
        if control == 'cost_control':
          
            AM = PM+IsM+IaM+VPM+VIsM+VIaM+WPM+WIsM+WIaM 
            nu_2 = numax*(max(C-Cc,0)/(max(C-Cc,0)+C0-Cc))*max(eta*Mc-AM,0)
        elif control == 'vac_control':
            
            if t<300:
                 nu_2 = numax*(max(C-Cc,0)/(max(C-Cc,0)+C0-Cc))
            else:
                
                V_0 = 0.75*n ; Vc = 0.4*n
                nu_2 = numax*(max(C-Cc,0)/(max(C-Cc,0)+C0-Cc))*(max(Vac-Vc,0)/(max(Vac-Vc,0)+V_0 - Vc))
        else:
            nu_2 = numax*(max(C-Cc,0)/(max(C-Cc,0)+C0-Cc))
            
         
        
        
        return nu_2*2.
           
  
    dy = np.zeros(83)
    dy[0] = -(beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*S - mu(t)*S + 0.5*nu(t)*S1 + (1-q2)*nu(t)*S2 -p*S+wI*(Rs+RsM) + wA*(Ra+RaM) +w*WS#+ wp*(VS+WS)
    dy[1] = -delta*(beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*S1 - 0.5*mu(t)*S1 +q1*mu(t)*S - 0.5*nu(t)*S1 + q2*nu(t)*S2 - p*S1+wI*Rs1+wA*Ra1 +w*WS1 #+ wp*(VS1+WS1)
    dy[2] = (1-q1)*mu(t)*S +0.5*mu(t)*S1 -nu(t)*S2 - p*S2+wI*Rs2+wA*Ra2 +w*WS2
    dy[3] = (beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*S - mu(t)*E + 0.5*nu(t)*E1 + (1-q2)*nu(t)*E2 -sigma*E - p*E + w*WE
    dy[4] = delta*(beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*S1 - 0.5*mu(t)*E1 + q1*mu(t)*E - 0.5*nu(t)*E1 + q2*nu(t)*E2 - sigma*E1 -p*E1 + w*WE1
    dy[5] = (1-q1)*mu(t)*E + 0.5*mu(t)*E1 -nu(t)*E2 - sigma*E2 - p*E2 + w*WE2
    dy[6] = sigma*E -mu(t)*P + 0.5*nu(t)*P1 + (1.-q2)*nu(t)*P2 - phi*P - rho0*P - p*P + w*WP
    dy[7] = sigma*E1 -0.5*mu(t)*P1 + q1*mu(t)*P - 0.5*nu(t)*P1 + q2*nu(t)*P2 - phi*P1 - rho0*P1 - p*P1 + w*WP1
    dy[8] = sigma*E2 + (1-q1)*mu(t)*P + 0.5*mu(t)*P1 -0.5*nu(t)*P2 - phi*P2 - rho0*P2 - p*P2 + w*WP2
    dy[9] = rho0*(P+P1+P2) - phi*PM + w*WPM
     
        
    dy[10] = q*phi*P - muI*Is - gamma*Is - rhoI*Is + w*WIs
    dy[11] = q*phi*P1 + qI*muI*Is - gamma*Is1 - rhoI*Is1 + w*WIs1
    dy[12] = q*phi*P2 + (1-qI)*muI*Is - gamma*Is2 - rhoI*Is2  + w*WIs2
    dy[13] = rhoI*(Is+Is1+Is2) + q*phi*PM - gamma*IsM + w*WIsM
    
    dy[14] = (1-q)*phi*P - mu(t)*Ia + 0.5*nu(t)*Ia1 + (1-q2)*nu(t)*Ia2 - gamma*Ia - rho0*Ia - p*Ia + w*WIa
    dy[15] = (1-q)*phi*P1 - 0.5*mu(t)*Ia1 + q1*mu(t)*Ia - 0.5*nu(t)*Ia1 + q2*nu(t)*Ia2 - gamma*Ia1 - rho0*Ia1 - p*Ia1  + w*WIa1
    dy[16] = (1-q)*phi*P2 + (1-q1)*mu(t)*Ia +0.5*mu(t)*Ia1 - nu(t)*Ia2 - gamma*Ia2 - rho0*Ia2 - p*Ia2  + w*WIa2
    dy[17] = rho0*(Ia+Ia1+Ia2) +(1-q)*phi*PM - gamma*IaM  + w*WIaM
    dy[18] = gamma*Is - p*Rs-wI*Rs-mu(t)*Rs +0.5*nu(t)*Rs1 + (1-q2)*nu(t)*Rs2 +w*WRs
    dy[19] = gamma*Is1 - p*Rs1-wI*Rs1- 0.5*mu(t)*Rs1 + q1*mu(t)*Rs -0.5*nu(t)*Rs1 + q2*nu(t)*Rs2 +w*WRs1
    dy[20] = gamma*Is2 - p*Rs2-wI*Rs2 +(1-q1)*mu(t)*Rs +0.5*mu(t)*Rs1 - nu(t)*Rs2 +w*WRs2
    dy[21] = gamma*IsM - p*RsM-wI*RsM +w*WRsM
    dy[22] = gamma*Ia - p*Ra-wA*Ra- mu(t)*Ra +0.5*nu(t)*Ra1 + (1-q2)*nu(t)*Ra2 + w*WRa
    dy[23] = gamma*Ia1 - p*Ra1-wA*Ra1 -0.5*mu(t)*Ra1 + q1*mu(t)*Ra - 0.5*nu(t)*Ra1 + q2*nu(t)*Ra2 + w*WRa1
    dy[24] = gamma*Ia2 - p*Ra2-wA*Ra2+(1-q1)*mu(t)*Ra +0.5*mu(t)*Ra1 - nu(t)*Ra2 + w*WRa2
    dy[25] = gamma*IaM - p*RaM-wA*RaM + w*WRaM
    
   
    dy[26] = rho0*(Ia+Ia1+Ia2+P+P1+P2 )+rhoI*(Is+Is1+Is2) + rhoV0*(VIa+WIa+VIa1+WIa1+VIa2+WIa2+VP+WP+VP1+WP1+VP2+WP2) + rhoVI*(VIs+WIs+VIs1+WIs1+VIs2+WIs2) 
   
    dy[27] = n*((S2+E2+VS2+VE2+WS2+WE2)+(1-delta)*(S1+E1+VS1+VE1+WE1+WS1) ) - muC*C
        
        
    dy[28] = -epsilon*(beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*VS - mumaxV*mu(t)*VS + 0.5*numaxV*nu(t)*VS1 + (1-q2)*numaxV*nu(t)*VS2 + qVS*p*S-wp*VS
    dy[29] = -epsilon*delta*(beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*VS1 - 0.5*mumaxV*mu(t)*VS1 + q1*mumaxV*mu(t)*VS - 0.5*numaxV*nu(t)*VS1 +\
               q2*numaxV*nu(t)*VS2 + qVS*p*S1-wp*VS1
    dy[30] = (1-q1)*mumaxV*mu(t)*VS + 0.5*mumaxV*mu(t)*VS1 - numaxV*nu(t)*VS2 + qVS*p*S2-wp*VS2
    dy[31] = epsilon*(beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*VS - mumaxV*mu(t)*VE + 0.5*numaxV*nu(t)*VE1 + (1-q2)*numaxV*nu(t)*VE2 -sigma*VE + qVE*p*E
    dy[32] = epsilon*delta*(beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*VS1 - 0.5*mumaxV*mu(t)*VE1 + q1*mumaxV*mu(t)*VE - 0.5*numaxV*nu(t)*VE1 + q2*numaxV*nu(t)*VE2 - sigma*VE1 + qVE*p*E1
    dy[33] = (1-q1)*mumaxV*mu(t)*VE + 0.5*mumaxV*mu(t)*VE1 - numaxV*nu(t)*VE2 - sigma*VE2 + qVE*p*E2
    dy[34] = sigma*VE -mumaxV*mu(t)*VP + 0.5*numaxV*nu(t)*VP1 + (1-q2)*numaxV*nu(t)*VP2 - phi*VP - rhoV0*VP + qVP*p*P
    dy[35] = sigma*VE1 -0.5*mumaxV*mu(t)*VP1 + q1*mumaxV*mu(t)*VP - 0.5*numaxV*nu(t)*VP1 + q2*numaxV*nu(t)*VP2 - phi*VP1 - rhoV0*VP1 + qVP*p*P1
    dy[36] = sigma*VE2 + (1-q1)*mumaxV*mu(t)*VP + 0.5*mumaxV*mu(t)*VP1 - numaxV*nu(t)*VP2 - phi*VP2 - rhoV0*VP2 + qVP*p*P2
    dy[37] = rhoV0*(VP+VP1+VP2) - phi*VPM
   
    dy[38] = qv*phi*VP - muI*VIs - gamma*VIs - rhoVI*VIs 
    dy[39] = qv*phi*VP1 + qI*muI*VIs - gamma*VIs1 - rhoVI*VIs1 
    dy[40] = qv*phi*VP2 + (1-qI)*muI*VIs - gamma*VIs2 - rhoVI*VIs2 
    dy[41] = rhoVI*(VIs+VIs1+VIs2) + qv*phi*VPM - gamma*VIsM 
    
        
    dy[42] = (1-qv)*phi*VP - mumaxV*mu(t)*VIa + 0.5*numaxV*nu(t)*VIa1 + (1-q2)*numaxV*nu(t)*VIa2 - gamma*VIa - rhoV0*VIa + qVIa*p*Ia
    dy[43] = (1-qv)*phi*VP1 - 0.5*mumaxV*mu(t)*VIa1 + q1*mumaxV*mu(t)*VIa - 0.5*numaxV*nu(t)*VIa1 + q2*numaxV*nu(t)*VIa2 - gamma*VIa1 - rhoV0*VIa1 + qVIa*p*Ia1
    dy[44] = (1-qv)*phi*VP2 + (1-q1)*mumaxV*mu(t)*VIa + 0.5*mumaxV*mu(t)*VIa1 - numaxV*nu(t)*VIa2 - gamma*VIa2 - rhoV0*VIa2 + qVIa*p*Ia2
    dy[45] = rhoV0*(VIa+VIa1+VIa2) +(1-qv)*phi*VPM - gamma*VIaM
    dy[46] = gamma*VIs + qVR*p*Rs-wV*VRs- mumax*mu(t)*VRs + 0.5*numax*nu(t)*VRs1 + (1-q2)*numax*nu(t)*VRs2
    dy[47] = gamma*VIs1 + qVR*p*Rs1-wV*VRs1 - 0.5*mumax*mu(t)*VRs1 + q1*mumax*mu(t)*VRs - 0.5*numax*nu(t)*VRs1 + q2*numax*nu(t)*VRs2
    dy[48] = gamma*VIs2 + qVR*p*Rs2-wV*VRs2 +(1-q1)*mumax*mu(t)*VRs + 0.5*mumax*mu(t)*VRs1 - numax*nu(t)*VRs2
    dy[49] = gamma*VIsM + qVR*p*RsM-wV*VRsM
    dy[50] = gamma*VIa + qVR*p*Ra-wV*VRa- mumax*mu(t)*VRa + 0.5*numax*nu(t)*VRa1 + (1-q2)*numax*nu(t)*VRa2
    dy[51] = gamma*VIa1 + qVR*p*Ra1-wV*VRa1- 0.5*mumax*mu(t)*VRa1 + q1*mumax*mu(t)*VRa - 0.5*numax*nu(t)*VRa1 + q2*numax*nu(t)*VRa2
    dy[52] = gamma*VIa2 + qVR*p*Ra2-wV*VRa2+(1-q1)*mumax*mu(t)*VRa + 0.5*mumax*mu(t)*VRa1 - numax*nu(t)*VRa2
    dy[53] = gamma*VIaM + qVR*p*RaM-wV*VRaM     
       
    dy[54] = -(beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*WS - mumaxV*mu(t)*WS + 0.5*numaxV*nu(t)*WS1 + (1-q2)*numaxV*nu(t)*WS2 +(1-qVS)*p*S+wV*(VRs+VRa+ VRsM+VRaM) +wI*(WRs+WRa+WRsM+WRaM) + wp*VS - w*WS
    dy[55] = -delta*(beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*WS1 - 0.5*mumaxV*mu(t)*WS1 + q1*mumaxV*mu(t)*WS - 0.5*numaxV*nu(t)*WS1 + \
                     q2*numaxV*nu(t)*WS2 + (1-qVS)*p*S1+wI*(WRs1+WRa1) + wV*(VRs1+VRa1)+ wp*VS1  - w*WS1
    dy[56] = (1-q1)*mumaxV*mu(t)*WS + 0.5*mumaxV*mu(t)*WS1 - numaxV*nu(t)*WS2 + (1-qVS)*p*S2+ wI*(WRs2+WRa2) + wV*(VRs2+VRa2) + wp*VS2  - w*WS2
    dy[57] = (beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*WS  - mumaxV*mu(t)*WE + 0.5*numaxV*nu(t)*WE1 + (1-q2)*numaxV*nu(t)*WE2 - sigma*WE\
                + (1-qVE)*p*E - w*WE
    dy[58] = delta*(beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*WS1- 0.5*mumaxV*mu(t)*WE1 + q1*mumaxV*mu(t)*WE - 0.5*numaxV*nu(t)*WE1 + \
                             q2*numaxV*nu(t)*WE2 - sigma*WE1 + (1-qVE)*p*E1 - w*WE1
    dy[59] = (1-q1)*mumaxV*mu(t)*WE + 0.5*mumaxV*mu(t)*WE1 - numaxV*nu(t)*WE2 - sigma*WE2 + (1-qVE)*p*E2 - w*WE2
    dy[60] = sigma*WE -mumaxV*mu(t)*WP + 0.5*numaxV*nu(t)*WP1 + (1-q2)*numaxV*nu(t)*WP2 - phi*WP - rhoV0*WP + (1-qVP)*p*P - w*WP
    dy[61] = sigma*WE1 - 0.5*mumaxV*mu(t)*WP1 + q1*mumaxV*mu(t)*WP - 0.5*numaxV*nu(t)*WP1 + q2*numaxV*nu(t)*WP2 - phi*WP1 - rhoV0*WP1 \
               + (1-qVP)*p*P1  - w*WP1
    dy[62] = sigma*WE2 + (1-q1)*mumaxV*mu(t)*WP + 0.5*mumaxV*mu(t)*WP1 - numaxV*nu(t)*WP2 - phi*WP2 - rhoV0*WP2 + (1-qVP)*p*P2  - w*WP2
    dy[63] = rhoV0*(WP+WP1+WP2) - phi*WPM  - w*WPM
    
    dy[64] = q*phi*WP - muI*WIs - gamma*WIs - rhoVI*WIs  - w*WIs
    dy[65] = q*phi*WP1 + qI*muI*WIs - gamma*WIs1 - rhoVI*WIs1- w*WIs1
    dy[66] = q*phi*WP2 + (1-qI)*muI*WIs - gamma*WIs2 - rhoVI*WIs2  - w*WIs2
    dy[67] = rhoVI*(WIs+WIs1+WIs2) + q*phi*WPM - gamma*WIsM  - w*WIsM

    dy[68] = (1-q)*phi*WP - mumaxV*mu(t)*WIa + 0.5*numaxV*nu(t)*WIa1 + (1-q2)*numaxV*nu(t)*WIa2 - gamma*WIa - rhoV0*WIa + (1-qVIa)*p*Ia - w*WIa
    dy[69] = (1-q)*phi*WP1 - 0.5*mumaxV*mu(t)*WIa1 + q1*mumaxV*mu(t)*WIa - 0.5*numaxV*nu(t)*WIa1 + q2*numaxV*nu(t)*WIa2 - gamma*WIa1 \
                  - rhoV0*WIa1 + (1-qVIa)*p*Ia1 - w*WIa1
    dy[70] = (1-q)*phi*WP2 + (1-q1)*mumaxV*mu(t)*WIa + 0.5*mumaxV*mu(t)*WIa1 - numaxV*nu(t)*WIa2 - gamma*WIa2 - rhoV0*WIa2 + (1-qVIa)*p*Ia2 - w*WIa2
    dy[71] = rhoV0*(WIa+WIa1+WIa2) +(1-q)*phi*WPM - gamma*WIaM - w*WIaM
    dy[72] = gamma*WIs + (1-qVR)*p*Rs-wI*WRs - mumax*mu(t)*WRs + 0.5*numax*nu(t)*WRs1 + (1-q2)*numax*nu(t)*WRs2 - w*WRs
    dy[73] = gamma*WIs1 + (1-qVR)*p*Rs1-wI*WRs1- 0.5*mumax*mu(t)*WRs1 + q1*mumax*mu(t)*WRs - 0.5*numax*nu(t)*WRs1 + q2*numax*nu(t)*WRs2 - w*WRs1
    dy[74] = gamma*WIs2 + (1-qVR)*p*Rs2-wI*WRs2+(1-q1)*mumax*mu(t)*WRs + 0.5*mumax*mu(t)*WRs1 - numax*nu(t)*WRs2 - w*WRs2
    dy[75] = gamma*WIsM + (1-qVR)*p*RsM-wI*WRsM - w*WRsM
    dy[76] = gamma*WIa + (1-qVR)*p*Ra-wI*WRa- mumax*mu(t)*WRa + 0.5*numax*nu(t)*WRa1 + (1-q2)*numax*nu(t)*WRa2 - w*WRa
    dy[77] = gamma*WIa1 + (1-qVR)*p*Ra1-wI*WRa1- 0.5*mumax*mu(t)*WRa1 + q1*mumax*mu(t)*WRa - numax/2*nu(t)*WRa1 + q2*numax*nu(t)*WRa2 - w*WRa1
    dy[78] = gamma*WIa2 + (1-qVR)*p*Ra2-wI*WRa2+(1-q1)*mumax*mu(t)*WRa + 0.5*mumax*mu(t)*WRa1 - numax*nu(t)*WRa2 - w*WRa2
    dy[79] = gamma*WIaM + (1-qVR)*p*RaM-wI*WRaM - w*WRaM
    dy[80] = p*(S+ S1+ S2+E+E1+E2+P+P1+P2+Ia+Ia1+Ia2+Rs+Rs1+Rs2+RsM+Ra+Ra1+Ra2+RaM) 
    dy[81] = delta*(beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*S1 + (beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*S \
             + epsilon*delta*(beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*VS1 +epsilon*(beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*VS \
              +delta*(beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*WS1 + (beta_t(t)*F_s(t)+beta_v(t)*F_sv(t))*WS
    dy[82] =  (rho0*(S+ S1+ S2+E+E1+E2+P+P1+P2+Ia+Ia1+Ia2+Rs+Rs1+Rs2+RsM+Ra+Ra1+Ra2+RaM) + rhoI*(Is + Is1+Is2)+ rho0*(VS+ VS1+ VS2+VE+VE1+VE2+VP+VP1+VP2+VIa+VIa1+VIa2+VRs+VRs1+VRs2+VRsM+VRa+VRa1+VRa2+VRaM) + rhoI*(VIs + VIs1+VIs2)+rho0*(WS+ WS1+ WS2+WE+WE1+WE2+WP+WP1+WP2+WIa+WIa1+WIa2+WRs+WRs1+WRs2+WRsM+WRa+WRa1+WRa2+WRaM)  +  rhoI*(WIs + WIs1+WIs2))
    
    
   
    return dy

sf = int(sys.argv[1])
epsilon = float(sys.argv[2])
wI = float(sys.argv[3] ); wA = float(sys.argv[4]) ; wp = float(sys.argv[5]) ; wV = float(sys.argv[5]) ; w = float(sys.argv[6])
print(sys.argv[1])
data_pos = np.loadtxt("Data_pos.txt", dtype= float)
data_tot = np.loadtxt("Data_tot.txt", dtype = float)
data_vac = np.loadtxt("Data_vac.txt", dtype = float)

N_crit = PARAMS['N_crit'] 
data_pos = (data_pos - data_pos[0])/N_crit
data_tot = (data_tot - data_tot[0])/N_crit
data_vac = data_vac/N_crit


time_windows=  [0,89,163,290,324,367,393,425,457,476,492,540,587,619, 643]

len_t = len(time_windows)

def fit_func(t,v1,v2,v3,v4): 
     
    sol = odeint(old_model, y0,t, args = (v1,v2,v3,v4)) 
    
    M = sol[:,26] 
    Am = sol[:,9] + sol[:,13]+ sol[:,17]  +sol[:,37] + sol[:,41] + sol[:,45] + sol[:,63]+sol[:,67] + sol[:,71]
    vac = sol[:,80]
    

    return np.concatenate((Am,M,vac))   


K_cv = [] ; M_cv = [] ; rho_vec = [] ; p_vec = []   

N = PARAMS['N'] ;  T = PARAMS['time_sim'] ; N0 = PARAMS['N0'] ; I0 = 0.0002*N ; Sx0 = (0./100.)*N
I0_f = I0/N0; Sx0_f = Sx0/N0;  S0 = (N - I0_f - Sx0_f)/N0
y0 = np.zeros(83)

y0[0] = S0 ; y0[10] = I0_f

print(len(y0))

v_old =  [ 0.008, 0.00263 ,0.009, 0.0]

for j in range(len_t-1): 

    data_old = np.concatenate((data_pos[time_windows[j]:time_windows[j+1]],data_tot[time_windows[j]:time_windows[j+1]],\
                                data_vac[time_windows[j]:time_windows[j+1]]))
    time_old = np.arange(time_windows[j], time_windows[j+1])
    vnew, kcov= curve_fit(fit_func, time_old,data_old , v_old,  bounds=(0, [0.5,0.5,0.0095, 0.5]) )
    v1 = vnew[0] ; v2 = vnew[1] ; v3 = vnew[2] ;  v4 = vnew[3]
    sol_new = odeint(old_model, y0,time_old, args = (v1,v2,v3,v4)) 
    

    #################################################################################
    K_cv.append(v1);  M_cv.append(v2);  rho_vec.append(v3); 
    p_vec.append(v4)
    
    len_t = len(time_old)
    
    v_old = vnew ; y0 = sol_new[:,][-1] 



K_c2 = np.array(K_cv) ; M_c2 = np.array(M_cv) ; rho_vec2 = np.array(rho_vec) 
p_vec2 = np.array(p_vec)

file1 = "data/Kc_%d.npy" %(sf)
file2 = "data/MC_%d.npy" %(sf)
file3 = "data/rhoc_%d.npy" %(sf)
file4 = "data/pc_%d.npy" %(sf)


np.save( file1, K_c2)  
np.save(file2, M_c2)  
np.save(file3, rho_vec2)  
np.save( file4,p_vec2)