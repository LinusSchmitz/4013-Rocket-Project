%% This is my code

close all
clear all
clc

% SI UNITS used for this code.  Conversions from Imperial to SI system made
% as appropriate

t(1) = 0; % [s] initial time
rb(1) = 0; % [m] initial burn grain displacement
T_vendor = csvread('Customer Reported Data.csv');  % vendor provided data
T_exp1 = csvread('LOG00042.csv');  % data from CU experimental static fire


%% INPUTS
cstar_eff = .9; % [-], cstar efficiency
t_step = 0.03; % [s] time step
P_atm = 101325.; % [Pa] ambient pressure
%P_atm = 83491.; % [Pa] ambient pressure
%a = .000005; % [-] burn rate coefficient
a = .000069917; % [-] burn rate coefficient
%a = .00015
%n = 0.5; % [-] burn rate exponent
n = 0.321; % [-] burn rate exponent
cstar(1) = 1500.; % [m/s] characteristic velocity
h_grain = 1.505; % [in] motor grain height
r_grain_i = 0.177/2; % [in] motor grain inner radius
r_grain_o = 0.908/2; % [in] motor grain outer radius
r_throat = 0.123/2; % [in] throat radius
r_exit = 0.231/2; % [in] exit radius
Mass = 0.025; % [kg]

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
V_chamber = 2*pi*r_grain_i^2*h_grain; % [m^3]
j = 1;
while rb < (r_grain_o - r_grain_i) && rb < h_grain % while there is unburned grain remaining
    [A_burn(j), V_burn(j+1), V_chamber(j+1)] = burn_geometry(r_grain_i,h_grain,rb,r_grain_o); % [m] burn area, burn cavity volume
    Pc(j) = ((a * rho_p * A_burn(j) * cstar(j)) / (A_throat)).^((1)/(1-n))/1e6; % [MPa] chamber pressure
    burn_rate(j) = a*(Pc(j)*1e6)^n; % [m/s] burn rate
    rb = rb + burn_rate(j) * t_step; % [m] updates burn displacement
    
    delta_Vol(j) = (V_chamber(j+1)-V_chamber(j))/t_step; % [m^3/s] rate of change in burn cavity volume 
    delta_Vol(j) = (V_burn(j+1)-V_burn(j))/t_step; % [m^3/s] rate of change in burn cavity volume 
    [T_predicted(j),cstar(j+1), m_dot(j),Tc(j)] = thrust_calc(P_atm, Pc(j), A_exit, rho_p, burn_rate(j), A_burn(j), AR_sup, delta_Vol(j));
    cstar(j+1) = cstar(j+1)*cstar_eff; % [m/s]
    if j == 1
        t(j) = t_step;
    else
        t(j) = t(j-1) + t_step;
    end
    j = j+1;
end




%% Plotting data to compare

figure 
hold on 
plot(T_vendor(:,1),T_vendor(:,2), "LineWidth", 3)
title("Thrust Curve Comparisons")
xlabel("Time (s)")
ylabel("Thrust (lbf)")  
xlim([0,1.2])
ylim([0,16])

time = T_exp1(:,2)-T_exp1(1,2);
force = T_exp1(:,1);
cforce = force - min(force);
cforce = cforce; %*4.44822162;
cforce(1:3144) = [];
time(1:3144) = [];
cforce(100:end) = [];
time(100:end) = [];
time = time-time(1);
plot(time/1000, cforce, "LineWidth", 3)

plot(t,T_predicted, 'LineWidth',3)

hold off
legend('Spec Sheet', 'Real Data','Solution Code');

%% Comparing values
%Isp = T_predicted./m_dot;
I_pred = trapz(t,T_predicted)* 4.44822162;
Isp_pred = I_pred/(Mass*9.81);
I_Vend = trapz(T_vendor(:,1),T_vendor(:,2))* 4.44822162;
Isp_Vend = I_Vend/(Mass*9.81);
I_exp = trapz(time/1000,cforce) *  4.44822162;
Isp_exp = I_exp/(Mass*9.81);

% C (effective exhaust velocity, not characteristic velocity) comparison: 
C_pred = Isp_pred * 9.81;
C_Vend = Isp_Vend * 9.81;
C_exp = Isp_exp * 9.81;


%% aero validation 
Pc_test = [3:0.05:4];
Pc_test_en = Pc_test * 145.038;

for i = 1:length(Pc_test_en)
    Pc_testen = Pc_test_en(i);
    [Test1(i),Test2(i)] = aerothermochemistry(Pc_testen,AR_sup);
end

    Mach1 = [Test1.Mach]';
    Mach2 = [Test2.Mach]';
    P01 = [Test1.P0]'/10;
    P02 = [Test2.P0]'/10;
    Cstar1 = [Test1.Cstar]';
    Cstar2 = [Test2.Cstar]';
    Alpha1 = [Test1.a]';
    Alpha2 = [Test2.a]';
    Rho_g1 = [Test1.rho]';
    Rho_g2 = [Test2.rho]';
    Pe1 = [Test1.P]'/10;
    Pe2 = [Test2.P]'/10;
    TestData = [Pc_test_en',Mach1,P01,Cstar1,Alpha1,Rho_g1,Pe1,Mach2,P02,Cstar2,Alpha2,Rho_g2,Pe2];