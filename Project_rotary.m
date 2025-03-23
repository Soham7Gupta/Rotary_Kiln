
% Calculate temperature at the internal wall (Tw)
Cg = 1.2;%Average specific heat of gas
Av=0.07;%Absorptivity
Cs = 0.2;%Average Specific heat of solid
Cv = 2.1;%Average specific heat of steam
Cw = 0.7;
es = 0.1;% solid emissivity
ew = 0.9;%Emissivity of internal wall
hgs = 0.01;%Heat transfer coefficient between solid and gas
ht = 4.3;
hw = 0.001;
Lcu = 1.848;
Les = 2.356;
Lli = 3.927;
ew0=0.25; % Emissivity of outer wall made of iron
Kw=0.205;% Thermal Conductivity of Iron
sigma = 5.67*10^(-8);
ht = 4.3; 
Hv = 2260;
Vs = 1;  % Assumed value, needs verification
kc = 6.005;  % Assumed value, needs verification
hgw=15;   % Assumed value, needs verification
eg=0.96; 
C3 = hgs * Lcu;
C4 = (sigma * Lcu * es) / (1 - (1 - es) * (1 - Av));
C5 = hgw * Lli;
C7 = hw * Les;
phi_sw = 1 / (1 - eg - (1 - ew) * (Lcu / Lli * (1 - es) + (1 - Lcu / Lli)));
C8 = sigma * Lcu * phi_sw * ew * es;
%Certain values known
Re=0.2;
Ri=0.1;
Ta=300;
ho=75;
Tg=750;
e=2.71;
Qa=5*10^(-3);
Tw0=800;
Tw=1200;
hi = hgw + sigma * eg * ew0 * (Tg^3 - Tw^3);
m=Ta/Tw0;
numerator=(Ta/(hi*Ri)+Tg*(1/((1+m+m^2+m^3)*ew0*sigma*Tw0^3+ho)*Re)-(log(Ri/Re))/Kw);
denominator=(1/(hi*Ri)+1/((1+m+m^2+m^3)*ew0*sigma*Tw0^3+ho)*Re-(log(Ri/Re))/Kw);
Tw = numerator / denominator;
hi = hgw + sigma * eg * ew0 * (Tg^3 - Tw^3);
%%



%For moisture Content
z_span=[0,3];
Qh0=3.33*10^(-4);
function dQhdz=odefun(z,Qh)
   ht=4.3;
   Tg=750;
   Ts=500;
   Qs=3.33*10^(-3);
   Hv=2260000;
   A=20;
   dQhdz=-ht * A * (Tg - Ts) * Qh / (Hv * (0.1 * Qs));
end
[zsol,Qh_sol]=ode45(@odefun,z_span,Qh0);
figure;
plot(zsol,Qh_sol)
xlabel("z:Length of Kiln(m)")
ylabel("Qh:Moisture flow rate(Kg/s)")
title("Variation of Qh along the length of kiln")
%For solid
z_span=[0,3]
Qs0=12;
function dQsdz=odefun1(z,Qs)
     kc=6.005;
     e=2.72;
     Ts=500;
     Qa=5*10^(-3);%kg/s
     ht=4.3;
     Av=0.07;
     Tg=750;
     Ts=500;
     Hv=2260;
     Qh=3.33*10^(-3);
     Vs=((Qs/3.14*(0.3)^2)/2000)*0.9;
     A=20;
     dQsdz = -kc * e^-(8033 / Ts) * Qs * Qa / ((Qs/3.14*(0.3)^2)/2000)*0.9; - ht * A * (Tg - Ts) * Qh / (Hv * (0.1 * Qs));
end
[zsol1,Qs_sol]=ode45(@odefun1,z_span,Qs0);
figure;
plot(zsol1,Qs_sol)

xlabel("z:Length of Kiln(m)")
ylabel("Qs:Solid flow rate(Kg/s)")
title("Variation of Qs along the length of kiln")
% Define constants for the whole system
z_span = [0, 4];     % Span for the length of the kiln (m)
Qh0 = 3.33e-4;       % Initial moisture content (kg)
Ts0 = 200;           % Initial solid temperature (K)
Tg0 = 600;           % Initial gas temperature (K)

% Parameters (These can be adjusted or kept as constants)
ht = 4.3;            % Heat transfer coefficient for moisture evaporation (W/m^2 K)
Hv = 2260e3;         % Latent heat of vaporization (J/kg)
A = 10;              % Heat transfer area (m^2)
Cs = 0.2e3;          % Specific heat capacity of the solid (J/kg K)
sigma = 5.67e-8;     % Stefan-Boltzmann constant (W/m^2 K^4)
Lcu = 1.848;         % Characteristic length for convection (m)
Les = 2.356;         % Characteristic length for radiation (m)
hgs = 10;            % Heat transfer coefficient between gas and solid (W/m^2 K)
ew = 0.9;            % Wall emissivity
es = 0.1;            % Solid emissivity
Av = 0.07;           % Surface area factor for radiation
Tw = 1400;           % Wall temperature (K)
Cv = 2.1e3;          % Specific heat of the vapor (J/kg K)
Cg = 1.2e3;          % Specific heat of the gas (J/kg K)
Qg0 = 11.65e-3;      % Initial heat content of the gas (J/kg)
eg = 0.1;            % Emissivity of the gas

% Other assumed parameters
Kw = 0.205; 
Qa = 5e-3;
deltaH = 131e3; 
Qs0 = 3.33e-3;       % Heat content of the solid
Vs = ((Qs0 / (3.14 * (0.3)^2)) / 2000) * 0.9;  % Volume of the solid in kiln

% ODE for Tg (gas temperature)
function dTgdz = odefunTg(z, Tg, Ts, Qh, ht, Hv, A, Cv, Qg0, hgs, sigma, es, Av, Lcu, Tw, eg, ew, Les)
    % Convective and radiative heat transfer terms
    C3 = hgs * Lcu;  
    C4 = (sigma * Lcu * es) / (1 - (1 - es) * (1 - Av));  
    C7 = (sigma * Les * ew) / (1 - (1 - ew) * (1 - Av));

    % Rate of change of gas temperature Tg
    dTgdz = (-C3 * (Tg - Ts) - C4 * (Tg^4 * eg - Ts^4 * Av) - C7 * (Tg^4 * eg - Tw^4 * Av) ...
             - ht * A * (Tg - Ts) * Qh / Hv) / (Qg0 * Cv);
end

% ODE for Ts (solid temperature)
function dTsdz = odefunTs(z, Ts, Tg, Qh, ht, A, Hv, Cs, hgs, sigma, es, Av, Lcu, Tw, eg, Kw, Qa, deltaH, Vs, Qs0, ew, Les)
    % Convective and radiative heat transfer terms
    C3 = hgs * Lcu;  
    C4 = (sigma * Lcu * es) / (1 - (1 - es) * (1 - Av));  
    Lli = 3.927;
    hw = 0.001;
    C7 = hw * Les;
    phi_sw = 1 / (1 - eg - (1 - ew) * (Lcu / Lli * (1 - es) + (1 - Lcu / Lli)));
    C8 = sigma * Lcu * phi_sw * ew * es;

    % Rate of change of solid temperature Ts
    dTsdz = (C3 * (Tg - Ts) + C4 * (Tg^4 * eg - Ts^4 * Av) + C7 * (Tw - Ts) + C8 * (Tw^4 * eg - Ts^4 * Av) ...
            - Kw * exp(-8033 / Ts) * (Qs0 * Qa * deltaH) / Vs) / (Qs0 * Cs);
end

% ODE for Qh (moisture content)
function dQhdz = odefunQh(z, Qh, Tg, Ts, ht, Hv, A)
    % Rate of change of moisture content Qh
    dQhdz = -ht * A * (Tg - Ts) * Qh / Hv;  
end

% Define the system of coupled ODEs for Tg, Ts, and Qh
function dydz = kiln_odes(z, y, ht, Hv, A, Cv, Qg0, hgs, sigma, es, Av, Lcu, Tw, eg, ew, Les, Cs, Kw, Qa, deltaH, Vs, Qs0)
    % y = [Tg, Ts, Qh]
    Tg = y(1);  % Gas temperature
    Ts = y(2);  % Solid temperature
    Qh = y(3);  % Moisture content

    % ODEs for Tg, Ts, and Qh
    dTgdz = odefunTg(z, Tg, Ts, Qh, ht, Hv, A, Cv, Qg0, hgs, sigma, es, Av, Lcu, Tw, eg, ew, Les);
    dTsdz = odefunTs(z, Ts, Tg, Qh, ht, A, Hv, Cs, hgs, sigma, es, Av, Lcu, Tw, eg, Kw, Qa, deltaH, Vs, Qs0, ew, Les);
    dQhdz = odefunQh(z, Qh, Tg, Ts, ht, Hv, A);

    % Return the derivatives as a column vector
    dydz = [dTgdz; dTsdz; dQhdz];
end

% Initial conditions: Tg0, Ts0, Qh0
initial_conditions = [Tg0, Ts0, Qh0];

% Solver options for higher accuracy (tighter tolerances)
options = odeset('RelTol',1e-6,'AbsTol',1e-8);

% Solve the coupled system using ode15s for stiff problems
[z_sol, y_sol] = ode15s(@(z, y) kiln_odes(z, y, ht, Hv, A, Cv, Qg0, hgs, sigma, es, Av, Lcu, Tw, eg, ew, Les, Cs, Kw, Qa, deltaH, Vs, Qs0), z_span, initial_conditions, options);

% Extract solutions
Tg_sol = y_sol(:, 1);  % Gas temperature
Ts_sol = y_sol(:, 2);  % Solid temperature
Qh_sol = y_sol(:, 3);  % Moisture content

% Plot the results for Tg, Ts, and Qh with higher resolution
figure;

plot(z_sol, Tg_sol, 'k', 'LineWidth', 2);
xlabel("z: Length of Kiln (m)");
ylabel("Tg: Gas Temperature (K)");
title("Gas Temperature (Tg) along the Length of the Kiln");
grid on;
hold on
z2=linspace(0,4,1000);
g=1000+(290*(1-exp(-z2/0.25)));
%plot(z2,g)
% Define the file path (replace with actual file name and your username)
s1 = '/Users/professional_Soham/Documents/s1.csv';

% Read the CSV file using readtable (for both numeric and non-numeric data)
data = readtable(s1);

% Assuming the first column is X data and the second column is Y data
x = data{:, 1};   % First column as X-axis
y = data{:, 2};   % Second column as Y-axis

% Plot the data
plot(x, y, '-r');
xlabel('z');  % Add X-axis label
ylabel('Tg');  % Add Y-axis label
title('Data');   % Add a title to the plot
grid on;                  % Turn on the grid for better readability
legend('Modelled data', 'Practical data','Location','southeast')


hold off


% Define constants for the whole system
z_span = [0, 4];     % Span for the length of the kiln (m)
Qh0 = 3.33e-4;       % Initial moisture content (kg)
Ts0 = 200;           % Initial solid temperature (K)

% Parameters (These can be adjusted or kept as constants)
ht = 4.3;            % Heat transfer coefficient for moisture evaporation (W/m^2 K)
Hv = 2260e3;         % Latent heat of vaporization (J/kg)
A = 10;             % Heat transfer area (m^2)
Cs = 0.2e3;          % Specific heat capacity of the solid (J/kg K)
sigma = 5.67e-8;     % Stefan-Boltzmann constant (W/m^2 K^4)
Lcu = 1.848;         % Characteristic length for convection (m)
Les = 2.356;         % Characteristic length for radiation (m)
hgs = 10;            % Heat transfer coefficient between gas and solid (W/m^2 K)
ew = 0.9;            % Wall emissivity
es = 0.1;            % Solid emissivity
Av = 0.07;           % Surface area factor for radiation
Tw = 1400;           % Wall temperature (K)
Cv = 2.1e3;          % Specific heat of the vapor (J/kg K)
Cg = 1.2e3;          % Specific heat of the gas (J/kg K)
Qg0 = 11.65e-3;      % Initial heat content of the gas (J/kg)
eg = 0.1;           % Emissivity of the gas

% Other assumed parameters
Kw = 0.205; 
Qa = 5e-3;
deltaH = 131e3; 
Qs0 = 3.33e-3;  % Heat content of the solid
Vs = ((Qs0 / (3.14 * (0.3)^2)) / 2000) * 0.9;  % Volume of the solid in kiln

% Polynomial for Tg(z)
Tg_poly = @(z) 1000+(290*(1-exp(-z/0.25)));

% ODE for Qh (moisture content)
function dQhdz = odefunQh22(z, Qh, ht, Hv, A, Tg_poly)
    % Calculate Tg from the polynomial
    Tg = Tg_poly(z);
    Ts = 700; % Keep Ts constant initially

    % Rate of change of moisture content
    dQhdz = -ht * A * (Tg - Ts) * Qh / Hv;  
end

% ODE for solid temperature (Ts)
function dTsdz = odefunTs22(z, Ts, Qh, Tg_poly, ht, A, Hv, Cs, hgs, sigma, es, Av, Lcu, Tw, eg, Kw, Qa, deltaH, Vs, Qs0, ew, Les)
    % Compute Tg from the polynomial
    Tg = Tg_poly(z);
    
    % Compute radiative and convective terms
    C3 = hgs * Lcu;  
    C4 = (sigma * Lcu * es) / (1 - (1 - es) * (1 - Av));  
    Lli = 3.927;
    hw = 0.001;
    C7 = hw * Les;
    phi_sw = 1 / (1 - eg - (1 - ew) * (Lcu / Lli * (1 - es) + (1 - Lcu / Lli)));
    C8 = sigma * Lcu * phi_sw * ew * es;

    % Rate of change of solid temperature Ts
    dTsdz = (C3 * (Tg - Ts) + C4 * (Tg^4 * eg - Ts^4 * Av) + C7 * (Tw - Ts) + C8 * (Tw^4 * eg - Ts^4 * Av) ...
            - Kw * exp(-8033 / Ts) * (Qs0 * Qa * deltaH) / Vs) / (Qs0 * Cs);
end

% Solve for Qh first using the updated ODE with variable Tg(z)
[z_sol, Qh_sol] = ode45(@(z, Qh) odefunQh22(z, Qh, ht, Hv, A, Tg_poly), z_span, Qh0);

% Solve for Ts using the updated ODE with variable Tg(z) and interpolated Qh values
[z_sol2, Ts_sol] = ode45(@(z, Ts) odefunTs22(z, Ts, interp1(z_sol, Qh_sol, z), Tg_poly, ht, A, Hv, Cs, hgs, sigma, es, Av, Lcu, Tw, eg, Kw, Qa, deltaH, Vs, Qs0, ew, Les), z_span, Ts0);
figure;
% Plot the results for Ts
plot(z_sol2, Ts_sol, 'k', 'LineWidth', 2);
xlabel("z: Length of Kiln (m)");
ylabel("Ts: Temperature of Solid (K)");
title("Variation of Solid Temperature (Ts) along the Length of Kiln");
hold on

% Define the file path (replace with actual file name and your username)
g2 = '/Users/professional_Soham/Documents/s2.csv';

% Read the CSV file using readtable (for both numeric and non-numeric data)
data = readtable(g2);

% Assuming the first column is X data and the second column is Y data
x = data{2:107, 1};   % First column as X-axis
y = data{2:107, 2};   % Second column as Y-axis

% Plot the data
plot(x, y, '-r');
xlabel('Z');  % Add X-axis label
ylabel('Ts');  % Add Y-axis label
title('Data');   % Add a title to the plot
grid on;                  % Turn on the grid for better readability
legend('Modelled data', 'Practical data','Location','southeast')
hold off