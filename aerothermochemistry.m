function[SM_results1, SM_results3] = aerothermochemistry(Pc_en, AR_sup)
    SM_inputs = CEAinput();                     %define object class to be passed into CEA code

    % Define Propellant Ingredients & initial condition 
    
    % For example
    SM_inputs.ox1      = 'NH4CLO4(I)';           %primary oxidizer
    SM_inputs.ox1T     = 536;                    %primary ox temp (K) 298

    SM_inputs.fu1      = 'C4H6,butadiene';       %primary fuel
    SM_inputs.fu1T     = 536;                    %primary fuel temp (K) 298

    SM_inputs.fu2      = 'AL(cr)';               %secondary fuel
    SM_inputs.fu2T     = 536;                    %secondary fuel temp (K) 298
      
    SM_inputs.ox1wt    = 84.85;                     %primary oxidizer weight (by mass) in total 80 82.5
    SM_inputs.fu1wt    = 5.15;                     %primary fuel weight (by mass) in total 10 7.5 
    SM_inputs.fu2wt    = 10;                     %secondary fuel weight (by mass) in total 10 10
    
    SM_inputs.Pc       = Pc_en;                     %chamber pressure (Psia)
    
    SM_inputs.supar    = AR_sup;                 % nozzle expansion ratio 
     
    %% Insert Algorithms or Table of Theoretical Aerothermochemistry Data for the chamber, nozzle throat and nozzle exit
    
%         SM_results.T = 2616; % [K] chamber temperature
%         SM_results.g = 1.235; % ratio of specific heats
%         SM_results.MolWt = 22.959; % molecular weight
%         SM_results.rho = 7.276; % [kg/m3] gas density
%         SM_results.c = 1500; % [m/s] cstar

    %%  Utilize Matlab code for CEA found online & Verified against CEA online. 
        SM_inputs.runCEA();   
    
    %% Return results
    % All units will be in si as the output and input because of the 'si'
    SM_results1 = SM_inputs.getCEAresults(1,'si');
    %SM_results2 = SM_inputs.getCEAresults(2,'si');
    SM_results3 = SM_inputs.getCEAresults(3,'si');           %1 for chamber conditions, 2 for throat, 3 for nozzle exit: requires arg for area ratio or exit pressure), 'en' for english units
    