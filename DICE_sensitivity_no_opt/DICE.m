% clc
% clear;
Initialset;
% Initial year 2010; time step 5 yr

global  S S0 K0 Tat0 Tlo0 Mat0 Mlo0 Mup0 pbacktime pai theta2
global  L A T gama fval aopt

% Initialization
S = zeros(20, T);
y0 = A(1)*(K0^gama)*(L(1)/1000)^(1-gama);
S0 = [K0, Tat0, Tlo0, Mat0, Mlo0, Mup0, y0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];

%myoptions = optimset('Display','iter','FunValCheck','on','algorithm','sqp','MaxFunEvals',20000);

fval = DICE_fun(aopt);