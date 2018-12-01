function SimpleSuicide_SD_qi(handles)
global a0 a1 a2 ac ar aa Ic Ia Ir b_suicide
global N Tmax
%----------------------------------------------------------------------------------------------
% SD model for suicide
%----------------------------------------------------------------------------------------------
% parameters
% a1 = base rate of recovery program from suicidal
% a2 = internal growth rate (set to 1)
% ar = recovery program 
% ac = community program 
% aa = awareness program 
% a0 = base rate of suicidal tendency from PD
% b_suicide = rate of commit suicide from suicidal
% pn0= initial proportional PD in non-suicidal stock


% flow Rate constants
% internal growth                                = a2;                 
% community intervention                         = ac*Ic;              
% Sn with PD to suicidal (by awareness)          = a0*(1 - aa*Ia);     
% Ss to Sr with PD (by recovery)                 = a1*(1-ar*Ir);       
% Ss to Sr without PD (by recovery)              = a1*((ar)*Ir);       
% Sr no PD to Sn no PD                           = 1;                   
% Sr PD to Sn PD                                 = 1;                   
% sproportional of Ss to commit suicide          = b_suicide;           
%----------------------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------------
% parameters initializaion
a0 = 0.5;
a1 = 0.3; 
a2 = 0.1;
ar = 0.4;
ac = 0.6;
aa = 0.5;


Ir = 0.50;
Ic = 0.51;
Ia = 0.52;
b_suicide = 0.00001;

pn0 = 0.3;


% get the value by changing sliders 
if nargin
    a0 = get(handles.a0,'Value'); handles.a0_txt.String = sprintf('a0 = %.2f (base rate of suicidal tendency from PD)',a0);
    a1 = get(handles.a1,'Value'); handles.a1_txt.String = sprintf('a1 = %.2f (base rate of recovery program from suicidal)',a1);
    a2 = get(handles.a2,'Value'); handles.a2_txt.String = sprintf('a2 = %.2f (internal growth rate)',a2);
    ar = get(handles.ar,'Value'); handles.ar_txt.String = sprintf('ar = %.2f (recovery program coefficient)',ar);
    ac = get(handles.ac,'Value'); handles.ac_txt.String = sprintf('ac = %.2f (community program coefficient)',ac);
    aa = get(handles.aa,'Value'); handles.aa_txt.String = sprintf('aa = %.2f (awareness raising coefficient)',aa);
    Ir = get(handles.Ir,'Value'); handles.Ir_txt.String = sprintf('Ir = %.2f (recovery program index)',Ir);
    Ic = get(handles.Ic,'Value'); handles.Ic_txt.String = sprintf('Ic = %.2f (community program index)',Ic);
    Ia = get(handles.Ia,'Value'); handles.Ia_txt.String = sprintf('Ia = %.2f (awareness raising index)',Ia);
    pn0= get(handles.Pn,'Value'); handles.Pn_txt.String = sprintf('Pn0 = %.2f (initial proportion of PD)' ,pn0);
    b_suicide = get(handles.b_suicide,'Value'); handles.b_suicide_txt.String = sprintf('b_suicide = %.2f (rate of commit suicide from suicidal)',b_suicide);
end

SnnPD = pn0;
W0    = [1-SnnPD SnnPD 0 0 0 0];       % initial condition: all people in non-suicidal stock

Tmax  = 20;                            % maximum time
N     = 600;
tspan = linspace(0, Tmax, N);          % unit: day

% run the ODE model
[SnPD, SnnPD, Ss, SrPD, SrnPD, Sdeath] = rk_Suicide_SD( W0);
% combine the stocks into one matrix, for plot
Y = [SnPD;SnnPD;Ss;SrPD;SrnPD;Sdeath];
TLTS = {'S_n with PD',...
    'S_n without PD',...
    'S_s',...
    'S_r with PD',...
    'S_r without PD',...
    'S_{death}'};            % titles for each subplots

if nargin
    h_axes = {handles.axes1,handles.axes2,handles.axes3,...
        handles.axes4,handles.axes5,handles.axes6};
else
    fig = figure(110); fig.Name = 'SD_model'; clf, h_axes = [];
end

for k=1:size(Y,1)
    if nargin
        axes(h_axes{k}),
    else
        subplot(2,3,k),
    end
    plot(tspan,Y(k,:),'LineWidth',2),title(sprintf(TLTS{k})),
    % set different y limit in subplots
    if k <5
        ylim([0,1])
        if k > 2   
            ylim([0,0.6])
        end
    end
    if k ==5
        if max(Y(k,:)) <= 0.1
        ylim([0,0.1])
        else
        ylim([0,0.3]) 
        end
    end
    xlim([0,Tmax])  
    xlabel('time (year)')
end


function [SnPD, SnnPD, Ss,SrPD, SrnPD, Sdeath] = rk_Suicide_SD(W)
%----------------------------------------------------------------------------------------------
% Runge-Kutta solver for the suicide ODEs
% W = Sn_PD Sn_nPD Ss Sr_Pd Sr_nPd
%----------------------------------------------------------------------------------------------
global Tmax N                                            % Tmax = maximum time (months), N = number of time points

n    = size(W,1);                                        % number of samples
Tmin = 0;                                                % minimum time
h    = (Tmax - Tmin) / N;                                % time increment
SnPD    = zeros(n,N);
SnnPD   = zeros(n,N); Ss = SnPD;
SrPD    = SnnPD;   SrnPD = SnnPD;  Sdeath =  SnnPD;      % initialize the vecotors
for i = 1:N                                              % for each time step
    t       = Tmin + h*i;
    SnPD(:,i)    = W(:,1);                        % stock of non suicidal with PD
    SnnPD(:,i)   = W(:,2);                        % stock of suicidal without PD
    Ss(:,i)      = W(:,3);                        % stock of suicidal
    SrPD(:,i)    = W(:,4);                        % stock of recovery with PD
    SrnPD(:,i)   = W(:,5);                        % stock of recovery without PD
    Sdeath(:,i)  = W(:,6);                        % stock of death
    k1 = fdW(W);                                    % integrate one step using Runge-Kutta
    k2 = fdW(W + k1 * h/2);
    k3 = fdW( W + k2 * h/2);
    k4 = fdW( W + k3 * h);
    W  = W + (k1 + 2*k2 + 2*k3 + k4) * h/6;
end


function dW = fdW(W)
%----------------------------------------------------------------------------------------------
% COPEWELL ODE:
%   dW/dt = fdW(t,W)
%   W     = Sn_PD Sn_nPD Ss Sr_Pd Sr_nPd Sdeath
%----------------------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------------
global a0 a1 a2 ac ar aa Ic Ia Ir b_suicide          % parameters for stocks into AT

SnPD    = W(:,1);                        % stock of non suicidal with PD
SnnPD   = W(:,2);                        % stock of suicidal without PD
Ss      = W(:,3);                        % stock of suicidal
SrPD    = W(:,4);                        % stock of recovery with PD
SrnPD   = W(:,5);                        % stock of recovery without PD
Sdeath  = W(:,6);                        % stock of death


% rate constant
nnPD_nPD_rate      = a2;                 % internal growth
nPD_nnPD_rate      = ac*Ic;              % community intervention
nPD_s_rate         = a0*(1 - aa*Ia);     % Sn with PD to suicidal (by awareness)
s_rPD_rate         = a1*(1-ar*Ir);       % Ss to Sr with PD (by recovery)
s_rnPD_rate        = a1*((ar)*Ir);       % Ss to Sr without PD (by recovery)
rnPD_nnPD_rate     = 1;                  % Sr no PD to Sn no PD
rPD_nPD_rate       = 1;                  % Sr PD to Sn PD
s_suicide_rate     = b_suicide;          % proportional of Ss to commit suicide             
death_birth_rate   = 1;
% flow
nnPD_nPD_flow      = nnPD_nPD_rate .* SnnPD;   % flow from non Sn,noPD to non Sn,PD
nPD_nnPD_flow      = nPD_nnPD_rate .* SnPD ;   % flow from non Sn,PD to non Sn,noPD
nPD_s_flow         = nPD_s_rate    .* SnPD;    % flow from non Sn,PD to Ss
s_rPD_flow         = s_rPD_rate    .* Ss;      % flow from non Ss,PD to Sr,PD
s_rnPD_flow        = s_rnPD_rate   .* Ss;      % flow from non Ss,PD to Sr,noPD
rnPD_nnPD_flow     = rnPD_nnPD_rate.* SrnPD;   % flow from non Sr,noPD to Sn,noPD
rPD_nPD_flow       = rPD_nPD_rate  .* SrPD;    % flow from non Sr,PD to Sn,PD
s_suicide_flow     = s_suicide_rate.* Ss;      % flow from non Ss to death
death_birth_flow   = death_birth_rate.* Sdeath;% flow from death to Sn,noPD



% Derivatives
% temporary Derivatives
dSnnPD = -nnPD_nPD_flow + nPD_nnPD_flow + rnPD_nnPD_flow + death_birth_flow;            % dSnnPD/dt
dSnPD  = -nPD_nnPD_flow - nPD_s_flow    + nnPD_nPD_flow  + rPD_nPD_flow;
dSs    = -s_rPD_flow    - s_rnPD_flow   - s_suicide_flow + nPD_s_flow;
dSrnPD = -rnPD_nnPD_flow + s_rnPD_flow;
dSrPD  = - rPD_nPD_flow  + s_rPD_flow ;
ddeath =  - death_birth_flow + s_suicide_flow ;

% check if the temporary derivatives are feasible (net flow out is less than the stock). 
% if not, than there is no outflow from the stock 
if dSnnPD + SnnPD <= 0
    nnPD_nPD_flow = 0;
end
if dSnPD + SnPD <= 0
    nPD_nnPD_flow = 0;
    nPD_s_flow = 0;
end
if dSs + Ss <= 0
    s_rPD_flow = 0;
    s_rnPD_flow = 0;
    s_suicide_flow = 0;
end
if SrPD + dSrPD <= 0
    rPD_nPD_flow = 0;
end
if SrnPD + dSrnPD <= 0
    rnPD_nnPD_flow = 0;
end

if Sdeath + ddeath <= 0
    death_birth_flow = 0;
end

% recalculate Derivatives 
dSnnPD = -nnPD_nPD_flow + nPD_nnPD_flow + rnPD_nnPD_flow + death_birth_flow;            % dSnnPD/dt
dSnPD  = -nPD_nnPD_flow - nPD_s_flow    + nnPD_nPD_flow  + rPD_nPD_flow;
dSs    = -s_rPD_flow    - s_rnPD_flow   - s_suicide_flow + nPD_s_flow;
dSrnPD = -rnPD_nnPD_flow + s_rnPD_flow;
dSrPD  = - rPD_nPD_flow  + s_rPD_flow ;
ddeath =  - death_birth_flow + s_suicide_flow ;


dW(:,1) = dSnPD;
dW(:,2) = dSnnPD;
dW(:,3) = dSs;
dW(:,4) = dSrPD;
dW(:,5) = dSrnPD;
dW(:,6) = ddeath;
