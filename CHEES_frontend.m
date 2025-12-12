
function [pars,state] = CHEES_frontend(time,input_C,d13C_input,d13C_CO2,d13C_methane,methane_oxygen,k_CH4_diff,k_CH4_from_ocean,k_CH4_oxiwatercolum)
% function [pars,state] = MBOX_frontend(time,input_C,d13C_input,d13C_CO2)
%%%%%% MBOX 4 box ocean model Ruoyuan Qiu 2024 modified from Zhao et al.,2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Define parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% set up global structures
global stepnumber
global pars
global workingstate

pars.time = time;
pars.input_C = input_C;
pars.d13C_input = d13C_input;
pars.d13C_CO2 = d13C_CO2;

pars.d13C_methane = d13C_methane;
pars.methane_oxygen = methane_oxygen;

pars.k_CH4_diff = k_CH4_diff; % 0.1 - 0.5 
pars.k_CH4_from_ocean = k_CH4_from_ocean; % 0.2 - 0,.5
pars.k_CH4_oxiwatercolum = k_CH4_oxiwatercolum; % 0.1 - 0.5

%%%%%% water reservoir sizes in m3 (m=margins, s=surface, h= hi-lat, d=deep)
pars.vol_p  = 2.6e15 ;  %%%% approx volume of all shelves and slope to depth 100m, area pecentage 5%
pars.vol_di = 5.4e15 ;  %%%% approx volume of all shelves and slope in depth 100-1000m, area pecentage 5%
pars.vol_s  = 2.75e16 ; %%%% approx volume of suface water to depth 100m, area pecentage 76.5%
pars.vol_h  = 1.22e16 ; %%%% approx volume of hi-lat to depth 250m, area pecentage 13.5%
pars.vol_d  = 1.35e18 ;
pars.vol_ocean = pars.vol_p + pars.vol_di + pars.vol_s + pars.vol_h + pars.vol_d ;

%%%% mixing coefficient (Sv)
pars.mixcoeff_dip = 30.28;
pars.mixcoeff_ds  = 46.33;
pars.mixcoeff_dh  = 54.9;

%%%%%% inorganic carbon reservoirs in moles C
pars.CO2_a_0  = 5e16 ;% default 5e16 
pars.DIC_p_0  = 5.2e15 ;
pars.DIC_di_0 = 1.08e16 ;
pars.DIC_s_0  = 5.37e16 ;
pars.DIC_h_0  = 2.71e16 ;
pars.DIC_d_0  = 3e18 ;
pars.ALK_p_0  = 5.2e15 ;
pars.ALK_di_0 = 1.08e16 ;
pars.ALK_s_0  = 5.37e16 ;
pars.ALK_h_0  = 2.71e16 ;
pars.ALK_d_0  = 3e18 ;

%%%%%% C isotope composition
pars.d13c_atm_0    = -7 ;
pars.d13c_DIC_p_0  = 0.1 ;
pars.d13c_DIC_di_0 = 0.1 ;
pars.d13c_DIC_s_0  = 0.1 ;
pars.d13c_DIC_h_0  = 0.1 ;
pars.d13c_DIC_d_0  = 0.1 ;

%%%%%% POC reservoirs in moles C
pars.POC_p_0  = 630e12 ;
pars.POC_di_0 = 250e12 ;
pars.POC_s_0  = 2329e12 ;
pars.POC_h_0  = 1084e12 ;
pars.POC_d_0  = 56000e12 ;


%%%%%% Dissolved phosphate in moles
pars.DP_p_0  = 1.82e12 ;
pars.DP_di_0 = 7.56e12 ;
pars.DP_s_0  = 0.55e12 ;
pars.DP_h_0  = 16.27e12 ;
pars.DP_d_0  = 2970e12 ;


pars.d13c_POC_p_0  = -26; % -26
pars.d13c_POC_di_0 = -26;
pars.d13c_POC_s_0  = -26 ;
pars.d13c_POC_h_0  = -26 ;
pars.d13c_POC_d_0  = -26 ;


%%%%%% initial amount of O2 in moles C
pars.O2_a_0  = 3.7e19 ; % 3.7e19
pars.O2_p_0  = 6.705e14 ;
pars.O2_di_0 = 8.964e14 ;
pars.O2_s_0  = 9.139e15 ;
pars.O2_h_0  = 4.02e15 ;
pars.O2_d_0  = 1.823e17 ;

pars.O2_conc_p_0 = pars.O2_p_0/pars.vol_p ;
pars.O2_conc_di_0 =pars.O2_di_0/pars.vol_di;
pars.O2_conc_s_0 =pars.O2_s_0/pars.vol_s;
pars.O2_conc_h_0 =pars.O2_h_0/pars.vol_h;
pars.O2_conc_d_0 =pars.O2_d_0/pars.vol_d;
pars.Sredu_part = 0.0860;

%%%%%% initial amount of FeIII in moles
pars.FeIII_p_0  = 1.56e9 ;
pars.FeIII_di_0 = 3.24e9 ;
pars.FeIII_s_0  = 9.625e9 ;
pars.FeIII_h_0  = 3.66e9 ;
pars.FeIII_d_0  = 810e9 ;

%%%%%% initial amount of sulfate in moles
%%%% 28 mM
% pars.SO4_p_0  = 7.28e16;
% pars.SO4_di_0 = 1.512e17;
% pars.SO4_s_0  = 7.7e17;
% pars.SO4_h_0  = 3.416e17;
% pars.SO4_d_0  = 3.78e19;

%%%% 12 mM
pars.SO4_p_0  = 3.12e16;
pars.SO4_di_0 = 5.4e15*12;
pars.SO4_s_0  = 2.75e16*12;
pars.SO4_h_0  = 1.22e16 *12;
pars.SO4_d_0  = 1.35e18 *12;

%%%%%% initial amount of FeII in moles
pars.FeII_p_0  = 0;
pars.FeII_di_0 = 0;
pars.FeII_s_0  = 0;
pars.FeII_h_0  = 0;
pars.FeII_d_0  = 0;

%%%%%% initial amount of H2S in moles
pars.H2S_p_0  = 0;
pars.H2S_di_0 = 0;
pars.H2S_s_0  = 0;
pars.H2S_h_0  = 0;
pars.H2S_d_0  = 0;

%%%%%% present day rates
pars.k_ccdeg = 11.1e12 * 2 ; % 8e12
pars.k_carbw = 12e12 ; % 12e12
pars.k_sfw   = 0 ;       %%% Seafloor weathering
pars.k_mccb  = 20e12;  % pars.k_carbw + pars.k_ccdeg - pars.k_sfw  20e12
pars.k_silw  = pars.k_mccb - pars.k_carbw ;  %硅酸盐风化速率
basfrac = 0.3 ;

pars.k_granw = pars.k_silw * (1 - basfrac) ;
pars.k_basw = pars.k_silw * basfrac ;

%%%%%% organic C cycle
pars.k_ocdeg = 1.25e12 * 2  ;
pars.k_locb = 2.5e12 ; % 4.5e12 ; 2.5e12 
pars.k_mocb = 7e12  ;  % 4.5e12 ;
pars.k_oxidw = pars.k_mocb + pars.k_locb - pars.k_ocdeg ;

%%%%%% present P, Fe, pyrite and sulfate weathering rate
pars.k_phosw = 1.2088e+11 ; % 0.0967e12 ; %  0.0967e12 用系数改变作为强迫 1.2088e+11
pars.k_FeIIIw = 2e9;
pars.k_pyritew = 1.85e12;
pars.k_sulfatew = 1.25e12;

%%%%%% present pyrite and sulfate burial rate
pars.k_pyriteb_p  = 1.4e12;
pars.k_pyriteb_di = 0.45e12;
pars.k_sulfateb   = 1.25e12;

%%%%%% Redfeild ratio
pars.Red_C_P = 106;
pars.Red_C_N = 106/16;
pars.Red_C_O = 106/138;
pars.Red_C_Fe = 106*2000;
pars.Red_Fe_P = pars.Red_C_P/pars.Red_C_Fe;

pars.BE_p  = 6/40;
pars.BE_di = 6/40;
pars.BE_d  = 1/10;


%%%%%% Monod constant; mol/m3
pars.KP = 0.1e-3;
pars.KFe = 0.1e-6;
pars.KmO2 = 10e-3;  %%%for sulfate reduction in water column
pars.KmO2s = 10e-3;  %%%for sulfate reduction in sediments 10e-3
pars.KmO2m = 1e-3;  %%%for methane oxidation in water column
pars.KmFeIII = 10;
pars.KmSO4 = 0.5; %%%0.5,1.6


%%%%%% reaction rate constant m3 mol-1 yr-1
pars.kpy = 0.3708/1e3*24*365.25;
pars.kironO = 1.4e5;
pars.ksulfO = 1.6e2; %%%1.6e2
pars.kSironR = 8;
pars.kSide = 4e3; %%%mol m-3 yr-1

%%%%%% Ksp
pars.Kspside = 10^(-8.4)*1e6;%%%mol2 m-6
pars.KspFeSaq = 10^(-5.08)*1e3;%%%mol m-3

pars.STFeSaq = 10^(-5.7)*1e3;%%%mol m-3

pars.FeIIIa_s = 1.443e9; %%%
pars.FeIIIa_h = 0;  %%%% 0.00795e9
pars.FeIIhydro = 13.5e9;

%%%%%% initial amount of CH4 in moles C
pars.CH4_a_0 =  8e13; %   
pars.f_CH4_a = 10e12; %  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Initialise   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% set maximum step size for solver
% options = odeset('maxstep',1e3) ;
options = odeset('maxstep',1e3,'NonNegative',[1:11 18:27 33:59]) ;

%%%% run beginning
% fprintf('Beginning run: \n')

%%%% set stepnumber to 1
stepnumber = 1 ;

%%%%%%% model timeframe in years (0 = present day)
pars.whenstart = -3e6 ;
pars.whenend = 10000 ; % 0.33e6 ;

%%%% model start state
pars.startstate(1) = pars.CO2_a_0 ;
pars.startstate(2) = pars.DIC_p_0 ;
pars.startstate(3) = pars.DIC_di_0 ;
pars.startstate(4) = pars.DIC_s_0 ;
pars.startstate(5) = pars.DIC_h_0 ;
pars.startstate(6) = pars.DIC_d_0 ;

pars.startstate(7) = pars.ALK_p_0 ;
pars.startstate(8) = pars.ALK_di_0 ;
pars.startstate(9) = pars.ALK_s_0 ;
pars.startstate(10) = pars.ALK_h_0 ;
pars.startstate(11) = pars.ALK_d_0 ;

pars.startstate(12) = pars.CO2_a_0 * pars.d13c_atm_0 ;
pars.startstate(13) = pars.DIC_p_0 * pars.d13c_DIC_p_0 ;
pars.startstate(14) = pars.DIC_di_0 * pars.d13c_DIC_di_0 ;
pars.startstate(15) = pars.DIC_s_0 * pars.d13c_DIC_s_0 ;
pars.startstate(16) = pars.DIC_h_0 * pars.d13c_DIC_h_0 ;
pars.startstate(17) = pars.DIC_d_0 * pars.d13c_DIC_d_0 ;

pars.startstate(18) = pars.POC_p_0;
pars.startstate(19) = pars.POC_di_0;
pars.startstate(20) = pars.POC_s_0;
pars.startstate(21) = pars.POC_h_0;
pars.startstate(22) = pars.POC_d_0;

pars.startstate(23) = pars.DP_p_0;
pars.startstate(24) = pars.DP_di_0;
pars.startstate(25) = pars.DP_s_0;
pars.startstate(26) = pars.DP_h_0;
pars.startstate(27) = pars.DP_d_0;

pars.startstate(28) = pars.d13c_POC_p_0*pars.POC_p_0;
pars.startstate(29) = pars.d13c_POC_di_0*pars.POC_di_0;
pars.startstate(30) = pars.d13c_POC_s_0*pars.POC_s_0;
pars.startstate(31) = pars.d13c_POC_h_0*pars.POC_h_0;
pars.startstate(32) = pars.d13c_POC_d_0*pars.POC_d_0;

pars.startstate(33) = pars.O2_a_0;
pars.startstate(34) = pars.O2_p_0;
pars.startstate(35) = pars.O2_di_0;
pars.startstate(36) = pars.O2_s_0;
pars.startstate(37) = pars.O2_h_0;
pars.startstate(38) = pars.O2_d_0;

pars.startstate(39) = pars.FeIII_p_0;
pars.startstate(40) = pars.FeIII_di_0;
pars.startstate(41) = pars.FeIII_s_0;
pars.startstate(42) = pars.FeIII_h_0;
pars.startstate(43) = pars.FeIII_d_0;

pars.startstate(44) = pars.SO4_p_0;
pars.startstate(45) = pars.SO4_di_0;
pars.startstate(46) = pars.SO4_s_0;
pars.startstate(47) = pars.SO4_h_0;
pars.startstate(48) = pars.SO4_d_0;

pars.startstate(49) = pars.FeII_p_0;
pars.startstate(50) = pars.FeII_di_0;
pars.startstate(51) = pars.FeII_s_0;
pars.startstate(52) = pars.FeII_h_0;
pars.startstate(53) = pars.FeII_d_0;

pars.startstate(54) = pars.H2S_p_0;
pars.startstate(55) = pars.H2S_di_0;
pars.startstate(56) = pars.H2S_s_0;
pars.startstate(57) = pars.H2S_h_0;
pars.startstate(58) = pars.H2S_d_0;

pars.startstate(59) = pars.CH4_a_0;



%%%%%%% run the system  
[rawoutput.T,rawoutput.Y] = ode15s(@CHEES_equations,[pars.whenstart pars.whenend],pars.startstate,options);

%%%%% start time counter
% tic



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Postprocessing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% takes 'workingstate' from model and turns into 'state' %%%%%%%

%%%% size of output 
pars.output_length = length(rawoutput.T) ;
%%%%%%%%%% model finished output to screen
% fprintf('Integration finished \t') ; fprintf('Total steps: %d \t' , stepnumber ) ; fprintf('Output steps: %d \n' , pars.output_length ) 
% toc

%%%%%%%%% print final model states using final state for each timepoint
%%%%%%%%% during integration
% fprintf('assembling state vectors... \t')
% tic

%%%% trecords is index of shared values between ode15s output T vector and
%%%% model recorded workingstate t vector
[sharedvals,trecords] = intersect(workingstate.time,rawoutput.T,'stable') ;

%%%%%% assemble output state vectors
field_names = fieldnames(workingstate) ;
for numfields = 1:length(field_names)
    eval([' state.' char( field_names(numfields) ) ' = workingstate.' char( field_names(numfields) ) '(trecords) ; '])
end

%%%%%% done message
% fprintf('Done: \n ')
% endtime = toc ;
% fprintf('time (s): %d \n', endtime )

end