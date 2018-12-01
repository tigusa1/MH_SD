function SimpleSuicide_SD_qi(handles)
global a0 a1 a2 ac ar aa Ic Ia Ir b_suicide
global nt Tmax
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
%----------------------------------------------------------------------------------------------
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
% parameters initializaion
a0 = 0.5;
a1 = 0.3;
a2 = 0.1;
ar = 0.4;
ac = 0.6;
aa = 0.5;

Ir = 0.30;
Ic = 0.30;
Ia = 0.30;
b_suicide = 0.01;



pn0 = 0.3;

% change of I at time tchange
tchange = 2;              % time of change (years)
Ic_new  = 0.9;            % change I
Ir_new  = Ir;
Ia_new  = Ia;

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
    Ir_new = get(handles.Ir_new,'Value'); handles.Ir_new_txt.String = sprintf('Ir_new = %.2f (new recovery program index)',Ir_new);
    Ic_new = get(handles.Ic_new,'Value'); handles.Ic_new_txt.String = sprintf('Ic_new = %.2f (new community program index)',Ic_new);
    Ia_new = get(handles.Ia_new,'Value'); handles.Ia_new_txt.String = sprintf('Ia_new = %.2f (new awareness raising index)',Ia_new);
    trecovery = get(handles.trecovery,'Value'); handles.trecovery_txt.String = sprintf('trecovery = %.2f ()',trecovery);
    tchange = get(handles.tchange,'Value'); handles.tchange_txt.String = sprintf('tchange = %.2f ()',tchange);
    b_suicide = get(handles.b_suicide,'Value');
    handles.b_suicide_txt.String = sprintf('b_suicide = %.2f (rate of commit suicide from suicidal)',b_suicide);
end



Sn    = pn0;
W0    = [Sn 1-Sn 0 0 0 0];             % initial condition: all people in non-suicidal stock

Tmax  = 5;                             % maximum time
nt    = 600;
tscale= 12;                            % scale the time (month)
tspan = linspace(0, Tmax, nt)';        % unit: year
trecovery = 1/2;                       % time for recovery program (year)

% run the ODE model
[Sn, Snn, Ss, Sr, Srn, Sdeath] = rk_Suicide_SD(W0,tscale,trecovery,tchange,Ir_new,Ic_new,Ia_new);

% combine the stocks into one matrix, for plot
Y = [Sn Snn Ss Sr Srn Sdeath];
TLTS = {'S_n with PD',...
    'S_n without PD',...
    'S_s',...
    'S_r with PD',...
    'S_r without PD',...
    'S_{death}'};                      % titles for each subplots

if nargin
    h_axes = {handles.axes1,handles.axes2,handles.axes3,...
        handles.axes4,handles.axes5,handles.axes6};
else
    fig = figure(110); fig.Name = 'SD_model'; clf, h_axes = [];
end

for k=1:size(Y,2)
    if nargin
        axes(h_axes{k}),
    else
        subplot(2,3,k),
    end
    plot(tspan,Y(:,k)*100,'LineWidth',2),title(sprintf(TLTS{k})),
    % set different y limit in subplots
    if k <5
        ylim([0,100])
        if k > 2
            ylim([0,60])
        end
    end
    if k==5
        if max(Y(:,k)) <= 0.1
            ylim([0,10])
        else
            ylim([0,30])
        end
    end
    xlim([0,Tmax])
    xlabel('time (year)')
end


function [Sn, Snn, Ss,Sr, Srn, Sdeath] = rk_Suicide_SD( W0, tscale, trecovery,tchange,Ir_new,Ic_new,Ia_new )
%----------------------------------------------------------------------------------------------
% Runge-Kutta solver for the suicide ODEs
% W = Sn_ Sn_n Ss Sr_ Sr_n
%----------------------------------------------------------------------------------------------
global Tmax nt Ir Ic Ia                         % Tmax = maximum time (months), nt = number of time points

m    = size(W0,2);                              % number of variables
Tmin = 0;                                       % minimum time
h    = (Tmax - Tmin) / nt * tscale;             % time increment
W    = nan(nt,m);
dW   = W;
irecovery = ceil(trecovery*tscale/h);           % time steps for trecovery
ichange   = ceil(tchange  *tscale/h);           % time step of tchange
do_change = true;                               % command to change

W(1,:)= W0;                                     % initialize W

for i = 1:nt-1                                  % for each time step
    t          = Tmin + h*i;
    if do_change && i>ichange
        do_change = false;                      % turn off command to change
        Ir        = Ir_new;
        Ic        = Ic_new;
        Ia        = Ia_new;       
    end
    Wi = W(i,:);                                % get W at time t
    k1 = fdW(Wi);                               % integrate one step using Runge-Kutta
    k2 = fdW(Wi + k1 * h/2);
    k3 = fdW(Wi + k2 * h/2);
    k4 = fdW(Wi + k3 * h);
    dW(i,:) = (k1 + 2*k2 + 2*k3 + k4) * h/6;
    Wi1     = Wi + dW(i,:);
    if i>irecovery
        dWir     = dW(i-irecovery,:);           % flow at time t-trecovery
        Wi1(1:2) = Wi1(1:2) + dWir(4:5);
        Wi1(4:5) = Wi1(4:5) - dWir(4:5);
        Wi1(1)   = Wi1(1)   + dWir(6  );
        Wi1(6)   = Wi1(6)   - dWir(6  );
    end
    W(i+1,:)= Wi1;
end

Sn    = W(:,1);                        % stock of non suicidal with PD
Snn   = W(:,2);                        % stock of non suicidal without PD
Ss    = W(:,3);                        % stock of suicidal
Sr    = W(:,4);                        % stock of recovery with PD
Srn   = W(:,5);                        % stock of recovery without PD
Sdeath= W(:,6);                        % stock of death


function dW = fdW(W)
%----------------------------------------------------------------------------------------------
%   dW/dt = fdW(t,W)
%   W     = Sn_ Sn_n Ss Sr_ Sr_n Sdeath
%----------------------------------------------------------------------------------------------
global a0 a1 a2 ac ar aa Ic Ia Ir b_suicide          % parameters for stocks into AT

Sn    = W(1);                        % stock of non suicidal with PD
Snn   = W(2);                        % stock of non suicidal without PD
Ss    = W(3);                        % stock of suicidal
Sr    = W(4);                        % stock of recovery with PD
Srn   = W(5);                        % stock of recovery without PD
Sdeath= W(6);                        % stock of death

% rate constant
nn_n_rate      = a2;                 % internal growth
n_nn_rate      = ac*Ic;              % community intervention
n_s_rate       = a0*(1 - aa*Ia);     % Sn with PD to suicidal (by awareness)
s_r_rate       = a1*(1-ar*Ir);       % Ss to Sr with PD (by recovery)
s_rn_rate      = a1*((ar)*Ir);       % Ss to Sr without PD (by recovery)
rn_nn_rate     = 0;                  % Sr no PD to Sn no PD
r_n_rate       = 0;                  % Sr PD to Sn PD
s_suicide_rate = b_suicide;          % proportional of Ss to commit suicide
death_birth_rate = 0;

% flow
nn_n_flow        = nn_n_rate * Snn;   % flow from non Sn,noPD to non Sn,PD
n_nn_flow        = n_nn_rate * Sn ;   % flow from non Sn,PD to non Sn,noPD
n_s_flow         = n_s_rate  * Sn;    % flow from non Sn,PD to Ss
s_r_flow         = s_r_rate  * Ss;    % flow from non Ss,PD to Sr,PD
s_rn_flow        = s_rn_rate * Ss;    % flow from non Ss,PD to Sr,noPD
rn_nn_flow       = rn_nn_rate* Srn;   % flow from non Sr,noPD to Sn,noPD
r_n_flow         = r_n_rate  * Sr;    % flow from non Sr,PD to Sn,PD
s_suicide_flow   = s_suicide_rate* Ss;      % flow from non Ss to death
death_birth_flow = death_birth_rate* Sdeath;% flow from death to Sn,noPD

% Derivatives
% temporary Derivatives
dSnn = -nn_n_flow + n_nn_flow + rn_nn_flow + death_birth_flow;            % dSnn/dt
dSn  =  nn_n_flow - n_nn_flow - n_s_flow    + r_n_flow;
dSs  = -s_r_flow  - s_rn_flow + n_s_flow    - s_suicide_flow;
dSrn =            + s_rn_flow - rn_nn_flow;
dSr  =  s_r_flow              - r_n_flow;
ddeath = -death_birth_flow                  + s_suicide_flow;

% check if the temporary derivatives are feasible (net flow out is less than the stock).
% if not, than there is no outflow from the stock
if dSnn + Snn <= 0
    nn_n_flow = 0;
end
if dSn + Sn <= 0
    n_nn_flow = 0;
    n_s_flow = 0;
end
if dSs + Ss <= 0
    s_r_flow = 0;
    s_rn_flow = 0;
    s_suicide_flow = 0;
end
if Sr + dSr <= 0
    r_n_flow = 0;
end
if Srn + dSrn <= 0
    rn_nn_flow = 0;
end

if Sdeath + ddeath <= 0
    death_birth_flow = 0;
end

% recalculate Derivatives
dSnn = -nn_n_flow + n_nn_flow + rn_nn_flow + death_birth_flow;            % dSnn/dt
dSn  = -n_nn_flow - n_s_flow  + nn_n_flow  + r_n_flow;
dSs  = -s_r_flow  - s_rn_flow - s_suicide_flow + n_s_flow;
dSrn = -rn_nn_flow + s_rn_flow;
dSr  = -r_n_flow   + s_r_flow;
ddeath =  -death_birth_flow + s_suicide_flow;

dW    = nan(1,6);
dW(1) = dSn;
dW(2) = dSnn;
dW(3) = dSs;
dW(4) = dSr;
dW(5) = dSrn;
dW(6) = ddeath;
