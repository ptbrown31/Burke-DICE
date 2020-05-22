%Algorithm for calculating the social utility of each state under ceratin
%action (a3)

function [ U ]  =  Utility( S4, a4, t4 )
%Social Utility for a given state St
%   state = St
%   action = a
%   time = t2
global optlrsav sai1 sai2 sai3 theta1 theta2 L elasmu pai a0 dam_exp_coeff

if t4==1
    a4=a0;
end

%   Damage cost function
%Damage = 1-(sai1 * S4(2) + sai2 * (S4(2) ^ sai3)); % DICE-2013R
%Damage = exp(-1*dam_exp_coeff*S4(2));

% Damage = 1; % no damage

Damage = -0.132.*log(S4(2))-0.0428 + 1;
if Damage > 1 || isfinite(Damage) == 0; Damage = 1; end


%   Abatement cost = AbatedEmission * Carbonprice / theta2
Abate = (pai(t4) ^ (1 - theta2)) * theta1(t4) * (a4 ^ theta2);
%Abate = 0;

%   Net output after damage and abatement
Q = (Damage - Abate) * S4(7);

%   Consumption
C = (1 - optlrsav) * Q;

%   Consumption per capita
c = C / L(t4) * 1000;

%   Utility per capita
% u = 1 + c ^ (1 - alpha) / (1 - alpha); This is used in DICE-2007
u = (c^(1-elasmu)-1)/(1-elasmu)-1;

%what if we are just maximizing consumption per capita
%u = c;

%   Social utility
U = u * L(t4);
%U = u; %just use per capita utility (not weighing the future more because of more people)

end

