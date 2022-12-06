close all
clear all
clc

t(1) = 0; % [s] initial time
rb(1) = 0.177*0.5*.0254; % [m] initial burn grain displacement
T_claimed = csvread('Thrust.csv');

%% INPUTS
cstar_eff = 1; % [-], cstar efficiency
t_step = 0.05; % [s] time step
P_atm = 101325; % [Pa] ambient pressure
a = 5e-6; % [-] burn rate coefficient, from RP HW11, Purdue
n = 0.5; % [-] burn rate exponent, from RP HW11, Purdue
cstar = 1500 * cstar_eff; % [m/s] characteristic velocity 
h_grain = 1.505; % [in] motor grain height
r_grain = 0.908/2; % [in] motor grain radius
r_throat = 0.123/2; % [in] throat radius
r_exit = 0.231/2; % [in] exit radius
Mass = 0.025; % [kg]



%% CONVERSION
h_grain = h_grain*0.0254; % [m] motor grain height
r_grain = r_grain*0.0254; % [m] motor grain radius
r_throat = r_throat*0.0254; % [m] throat radius
r_exit = r_exit*0.0254; % [m] exit radius

%% QUANTITY CALCULATIONS
Vol = h_grain*r_grain^2*pi(); % [m^3]
rho_p = Mass/Vol; % [kg/m^3]
A_throat = pi()*(r_throat)^2; % [m^2]
A_exit = pi()*(r_exit)^2; % [m^2]
AR_sup = A_exit/A_throat; % supersonic area ratio

V_burn = 0; % [m^3]
j = 1;
while rb < r_grain && rb < h_grain % while there is grain remaining
    [A_burn(j), V_burn(j+1)] = burn_geometry2(r_grain,h_grain,rb); % [m] burn area, burn cavity volume
    Pc(j) = ((a * rho_p * A_burn(j) * cstar) / (A_throat)).^((1)/(1-n))/1e6; % [MPa] chamber pressure
    burn_rate(j) = a*(Pc(j)*10^6)^n; % [m/s] burn rate
    rb = rb + burn_rate(j) * t_step; % [m] updates burn displacement

    delta_Vol = (V_burn(j+1) - V_burn(j))/t_step; % [m^3/s] rate of change in burn cavity volume 
    [T_predicted(j),cstar,m_dot(j)] = thrust_calc(P_atm, Pc(j), A_exit, rho_p, burn_rate(j), A_burn(j), AR_sup, delta_Vol);
    cstar = cstar*cstar_eff; % [m/s]
    if j == 1
        t(j) = t_step;
    else
        t(j) = t(j-1) + t_step;
    end
    j = j+1
end

Isp = T_predicted./m_dot;

figure
hold on
plot(t,T_predicted)
plot(T_claimed(:,1),T_claimed(:,2))
title('Thrust Curve')
xlabel('Time (s)')
ylabel('Thrust (lbf)')

% trapz(t,T_predicted)
trapz(T_claimed(:,1),T_claimed(:,2))


fname = 'newtest.csv';
data = csvread(fname);
time = data(:,1) - data(1,1);
force = data(:,2);
cforce = force - min(force);
cforce = cforce * 3.3/6*25/255;
cforce(1:2900) = [];
time(1:2900) = [];
time = time-time(1);
plot(time/1000, cforce)
legend('Spec Sheet', 'Solution Code', 'Real Data');

% figure
% hold on
% plot(t,A_burn)