%% This is my code, their values 

close all
clear all
clc

% SI UNITS used for this code.  Conversions from Imperial to SI system made
% as appropriate

t(1) = 0; % [s] initial time
rb(1) = 0; % [m] initial burn grain displacement
T_vendor = csvread('Thrust.csv');  % vendor provided data
T_exp = csvread('Test.csv');  % data from CU experimental static fire

%% INPUTS
cstar_eff = 0.9; % [-], cstar efficiency
t_step = 0.05; % [s] time step
P_atm = 101325.; % [Pa] ambient pressure
a = .000004; % [-] burn rate coefficient
n = 0.48; % [-] burn rate exponent
cstar = 1500.; % [m/s] characteristic velocity
h_grain = 1.50; % [in] motor grain height
r_grain_i = 0.215/2; % [in] motor grain inner radius
r_grain_o = 0.93/2; % [in] motor grain outer radius
r_throat = 0.09/2; % [in] throat radius
r_exit = 0.2/2; % [in] exit radius
Mass = 0.0255; % [kg]

%% CONVERSIONS
h_grain = h_grain*0.0254; % [m] motor grain height
r_grain_i = r_grain_i*0.0254; % [m] motor grain inner radius
r_grain_o = r_grain_o*0.0254; % [m] motor grain outer radius
r_throat = r_throat*0.0254; % [m] throat radius
r_exit = r_exit*0.0254; % [m] exit radius

%% QUANTITY CALCULATIONS
Vol = h_grain*(r_grain_o^2-r_grain_i^2)*pi(); % [m^3]
rho_p = Mass/Vol; % [kg/m^3]
A_throat = pi()*(r_throat)^2; % [m^2]
A_exit = pi()*(r_exit)^2; % [m^2]
AR_sup = A_exit/A_throat; % supersonic area ratio

V_burn = 0; % [m^3]
%V_chamber = 0; % [m^3]
j = 1;
while rb < (r_grain_o - r_grain_i) && rb < h_grain % while there is unburned grain remaining
    [A_burn(j), V_burn(j), V_chamber(j)] = burn_geometry(r_grain_i,h_grain,rb); % [m] burn area, burn cavity volume
    Pc(j) = ((a * rho_p * A_burn(j) * cstar) / (A_throat)).^((1)/(1-n)); % [Pa] chamber pressure
    burn_rate(j) = a*(Pc(j))^n; % [m/s] burn rate
    rb = rb + burn_rate(j) * t_step; % [m] updates burn displacement
    
    % delta_Vol = ; % [m^3/s] rate of change in burn cavity volume 
    [T_predicted(j),cstar, m_dot(j)] = thrust_calc(P_atm, Pc(j), A_exit, rho_p, burn_rate(j), A_burn(j), AR_sup); %, delta_Vol);
    cstar = cstar*cstar_eff; % [m/s]
    if j == 1
        t(j) = t_step;
    else
        t(j) = t(j-1) + t_step;
    end
    j = j+1
end

Isp = T_predicted./m_dot;
