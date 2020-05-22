global T tstep gama sai1 sai2 sai3 deltak optlrsav pb pd theta1 theta2 pbacktime prstp
global Sig0 Sigg0 Siggd Sigg Sig E0 Eland Mat0 Mlo0 Mup0
global b11 b12 b21 b22 b23 b32 b33 F2010 F2100 Tat0 Tlo0 
global etha1 etha3 etha4 deltarf Fex L0 Lg0 LA L Ag Ag0 Agd A A0 pai ro
global R elasmu K0 lambda a0 deltaT cpricebase actbase alpha

% OUTPUT AND CAPITAL ACCUMULATION
%   T = Total Time horizon
T = 65;
time = 1:T;
tstep = 5; %  (5 years per period)
fosslim = 6000;      %  Maximum cumulative extraction fossil fuels (GtC) /6000/

% OUTPUT AND CAPITAL ACCUMULATION
%   capital share
gama = 0.3;
%   damage coefficient on temperature
sai1 = 0.0;
%   damage coefficient on temperature squared
sai2 = 0.0022; % DICE-2013R; 0.0028388 for DICE-2007
%   Exponent on damages
sai3 = 2.0;

%dam_exp_coeff = 0.1;

%   rate of depreciation (percent per year)
deltak = 0.1; %default = 0.1
%   Elasticity of marginal utility of consumption
%elasmu=1.45;
%elasmu=1.001;
%   Initial rate of social time preference per year
%prstp=0.015; %0.015;
%prstp=0.0; %0.015;

%   Optimal long-run savings rate used for transversality
optlrsav = (deltak + 0.004)/(deltak + 0.004*elasmu + prstp)*gama;
%   Initial capital stock ($ trillion 2005 USD)
K0 = 135; % DICE-2013R; 137 for DICE-2007

% EMISSIONS

% Industrial emissions 2010 (GtCO2 per year)
e0=35.85;
% Initial emissions control rate for base case 2010
miu0=0.03;
% Initial world gross output (trill 2005 USD)
q0=63.69;
%    Initial Sig
Sig0 = e0/(q0*(1-miu0)); % DICE-2013R
%    Initial growth rate of Sig (percent per decade)
%Sigg0 = -0.01; % DICE-2013R
%    Rate of decrease in the growth rate of Sig (percent per year)
%Siggd = -0.001; % DICE-2013R
%    Growth rate of Sig (percent per decade)
Sigg(1) = Sigg0;
%    Sig (industrial CO2 emissions/output -- MTC/$1000)
Sig = zeros(1,T);
Sig(1) = Sig0;
for i = 2:T
    Sigg(i) = Sigg(i-1)*((1+Siggd)^tstep);
    Sig(i) = Sig(i-1)*exp(Sigg(i)*tstep);
end
%override carbon intesntiy of output constant
%  Sig(:) = Sig0;
%  Sigg(:) = Sigg0;

%    Initial carbon emissions from land use change (GTCO2 per year)
E0 = 2.6; % 
%    Decline rate of land emissions (per period)
deland=0.2;
%    Carbon emissions from land use change (GTCO2 per year)
Eland = zeros(1,T);
for i=1:T % DICE-2013R
    Eland(i) = E0*(1-deland)^(time(i)-1);
end

%   Abatement cost function coefficient
theta1 = zeros(1,T+5); % theta1 = (pb*Sig/theta2).*((pr-1+exp(-pd*(time-1)))/pr); % DICE-2007
%   cost of backstop 2005 $ per t C in 2010



pb = 485; %3013: 344, 2016 550 (but in 2010 USD so ~485 in 2005 USD)
%   initial cost decline backstop cost per period
%pd = 0.025; % DICE-2013R, modified for a 5-year period; 0.05 for DICE-2007
%   Exponent of control cost function
%theta2 = 2.8;

%   Backstop price
pbacktime = zeros(1,T);
for i=1:T % DICE-2013R
    pbacktime(i) = pb*(1-pd)^(time(i)-1);
    theta1(i) = pbacktime(i)* (Sig(i) )/theta2/1000;
end
%override backstop price to make it constant
%  pbacktime(:) = pb;
%  theta1(:) = theta1(1);

% CONCENTRATIONS

%    Initial atmospheric concentration of CO2 (GTC) in 2010
Mat0 = 830.4;
%    Initial concentration of CO2 in biosphere/shallow oceans (GTC) in 2010
Mlo0 = 1527; % DICE-2013R; 1255 for DICE-2007 
%    Initial concentration of CO2 in deep oceans (GTC) in 2010
Mup0 = 10010; % DICE-2013R; 18365 for DICE-2007 
%    Carbon cycle transition coefficients (percent per decade)
%       atmosphere to atmosphere (b11)
b11 = 91.2; % DICE-2013R; 81.0712 for DICE-2007
%       biosphere/shallow oceans to atmosphere (b21)
b21 = 3.83; % DICE-2013R; 8.8 for DICE-2007
%        atmosphere to biosphere/shallow oceans (b12)
b12 = 8.80; % DICE-2013R; 18.9288 for DICE-2007
%       biosphere/shallow oceans to biosphere/shallow oceans (b22)
b22 = 95.92; % DICE-2013R; 85.2787 for DICE-2007
%        deep oceans to biosphere/shallow oceans (b32)
b32 = 0.03375; % DICE-2013R; 0.3119 for DICE-2007
%       biosphere/shallow oceans to deep oceans (b23)
b23 = 0.25; % DICE-2013R; 5 for DICE-2007
%       deep oceans to deep oceans (b33)
b33 = 99.97; % DICE-2013R; 99.6881 for DICE-2007

% TEMPERATURE

%   2010 forcings other ghg
F2010 = 0.25;
%   2100 forcings other ghg
F2100 = 0.70;
%   Initial atmospheric temperature (deg. C above 1900)
Tat0 = 0.8;
%   Initial temperature of deep oceans (deg. C above 1900)
Tlo0 = 0.0068;
%   Speed of adjustment parameter for atmospheric temperature
etha1 = 0.098; % DICE-2013R; 0.22 for DICE-2007
%   FCO22x Forcings of equilibrium CO2 doubling (Wm-2)
deltarf = 3.8;
%   Coefficient of heat loss from atmosphere to oceans
etha3 = 0.088; % DICE-2013R; 0.3 for DICE-2007
%   Coefficient of heat gain by deep oceans
etha4 = 0.025; % DICE-2013R; 0.05 for DICE-2007
%   Exogenous forcing (Watts per square meter)
Fex = zeros(1,T);
for i=1:T % DICE-2013R
    if i<19
        Fex(i) = F2010+ (1/18)*(F2100-F2010)*(time(i)-1);
    else
        Fex(i) = F2100;
    end
end

% POPULATION

%   Initial population (millions) 2010
L0 = 6838;
%   Growth rate to calibrate to 2050 pop projection
Lg0 = 0.134;
%   Asymptotic population
%LA = 10500;
%   Population (millions)
L = zeros(1,T);
L(1) = L0;
for i=1:(T-1) % DICE-2013R
    L(i+1) = L(i)*(LA/L(i))^Lg0;
end
%override to make constant
% L(:) = L(1);

% PRODUCTIVITY

%   Initial level of total factor productivity
A0 = 3.8; % DICE-2013R; 27.22 for DICE-2007
%   Initial growth rate for TFP per 5 years
Ag0 = 0.079; % DICE-2013R
%   Decline rate of TFP per 5 years
Agd = 0.006; % DICE-2013R
%   Rate of growth of productivity (percent per decade)
Ag = zeros(1,T);
%   Total factor productivity
A = zeros(1,T);
A(1) = A0;
for i = 1:(T-1)
    Ag(i)=Ag0*exp(-Agd*tstep*(time(i)-1));
    A(i+1) = A(i)/(1-Ag(i));
end
%override to make constant
%  Ag(:) = Ag(1);
%  A(:) = A(1);

% PARTICIPATION

%   PARTFRACT (participation rate)
pai = ones(1,T+5);

% WELFARE

%   Social rate of time preference (percent per year)
ro = prstp*100*ones(1,T);
%   Social time preference factor
R = zeros(1,T);
R(1) = 1;
for i = 2:T
    R(i) = R(i-1)/(1+ro(i-1)/100)^tstep;
end
%alpha = 2;
alpha = elasmu;
lambda = 1/(1+ro(1)/100)^tstep;

% Initialization

% Initial Action
a0 = 0.039;
%deltaT = 2.9; % Equilibrium temp impact (oC per doubling CO2)

% ** Abatement cost
% Period before which no emissions controls base
tnopol = 45;
% Initial base carbon price (2005$ per tCO2)
cprice0 = 1.0; 
% Growth rate of base carbon price per year
gcprice = 0.02; 
% *Base Case    Control Rate
actbase = zeros(1,T);
% *Base Case      Carbon Price
cpricebase = zeros(1,T);
for i = 1:T
    cpricebase(i) = cprice0*(1+gcprice)^(tstep*(time(i)-1));
    actbase(i) = pai(i) * exp( (log(cpricebase(i)/pbacktime(i))) / (theta2-1) );
end


