function [ThSM_en, Cstar, m_dot,Tc] = thrust_calc(Pa, Pc, Ae, rho_p, burn_rate, A_burn, AR_sup,deltaVol)
    Pc_en = Pc * 145.038; % [psia] chamber pressure

    % Pa in units of anything (just has to be consistent) 
    % Pc in units of Pa right now (needs to be psia for panning into
    % aerothermochem)
    % Ae exit area in units of [m^2] 
    % rho_p propellant density [kg/m^3]
    % burn_rate in units of [m/s]
    % A_burn in units of [m^2]
    % AR_sup nozel expansion ratio (unitless) 

    %% Aerothermochemistry Data
    ERROR = 0;
    try % tests if there is any output
        % OUTPUT1 GIVES VALUES IN THE CHAMBER
        % OUTPUT2 GIVES VALUES AT THE NOZZLE THROAT
        % OUTPUT3 GIVES VALUES AT NOZZLE EXIT
        [Output1, Output3] = aerothermochemistry(Pc_en, AR_sup); % Make sure all inputs are si
    catch
      ERROR = 1;
    end
    
    if ERROR == 1 % sets Mach and alpha to zero if output DNE
       alpha = 0;
       Mach = 0;
       Pe = Pa;
       Cstar = 0;
       rho_g = 0;
    else
        Tc = Output1.T; % [K] chamber temperature
        %gamma_c = Output1.GAMMA; % ratio of specific heats
        %MolWt_c = Output1.MolWt; % molecular weight
        rho_g = Output1.rho; % [kg/m3] gas density
        Cstar = Output1.Cstar; % [m/s] c* velocity

        Pe = Output3.P *1e5; % [Pa] exit pressure
        Mach = Output3.Mach; % exit Mach number
        alpha = Output3.a; % [m/s] exit speed of sound
        
    end

    %% MASS FLOW CALCULATION
    %m_dot = 0.000005; % [kg/s] propellant mass flow rate
    m_dot = (rho_p-rho_g)*A_burn*burn_rate;% - rho_g*deltaVol; % [kg/s] propellant mass flow rate

    %% THRUST CALCULATION
    %ThSM = 20; % [N] thrust SI
    ThSM = m_dot*(Mach*alpha)+(Pe - Pa)*Ae; % [N] thrust SI
    ThSM_en = ThSM * 0.224809; % [lbf] imperial thrust to match curve data