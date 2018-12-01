function SimpleSuicide_SD_qi(handles)
global a0  ac ar aa b0 Ic Ia Ir c
global N Tmax
%----------------------------------------------------------------------------------------------
% SD model for suicide
%----------------------------------------------------------------------------------------------
% a0 = internal growth rate (set to 1)
% ar = recovery program ES
% ac = community program ES
% aa = awareness program ES
% b0 = base rate of suicidal tendency from PD
% pn = proportion of non-suicidal that have PD
% pn0= baseline proportion
%----------------------------------------------------------------------------------------------
% CHANGE IN PD
% internal:  +a0*Sn*(pn0 - pn)
% community: -ac*Ic*Sn*pn
% recovery:  -ar*Ir*Sr
%
% STOCKS
% non-suicidal: Sn
% suicidal:     Ss = b*Sn*pn
% recovery:     Sr = Ss
%
% FLOW RATE
% suicidal:     b  = b0*(1 - aa*Ia)
%
% EQUILIBRIUM
% a0*(pn0 - pn) = ac*Ic*pn + ar*Ir*(b*pn)
%           pn  = a0*pn0/[a0 + ac*Ic + ar*Ir*b] = pn0 (if no I)
%
% pn*Sn = pn*(N - Ss - Sr) = pn*(N - 2*Ss) = pn*(N - 2*b*Sn*pn)
%    Sn = N - 2*b*Sn*pn
%    Sn = N / (1 + 2*b*pn)
%----------------------------------------------------------------------------------------------
a0 = 1;
ar = 0.4;
ac = 0.6;
aa = 0.5;
b0 = 0.1;

Ir = 0.50;
Ic = 0.51;
Ia = 0.52;
c = 0.00001;

pn0 = 0.3;

% number of time steps

tspan = linspace(0, Tmax, N); % unit: month

if nargin
    ar = get(handles.ar,'Value'); handles.ar_txt.String = sprintf('ar = %.2f',ar);
    ac = get(handles.ac,'Value'); handles.ac_txt.String = sprintf('ac = %.2f',ac);
    aa = get(handles.aa,'Value'); handles.aa_txt.String = sprintf('aa = %.2f',aa);
    b0 = get(handles.b0,'Value'); handles.b0_txt.String = sprintf('b0 = %.2f',b0);
    Ir = get(handles.Ir,'Value'); handles.Ir_txt.String = sprintf('Ir = %.2f',Ir);
    Ic = get(handles.Ic,'Value'); handles.Ic_txt.String = sprintf('Ic = %.2f',Ic);
    Ia = get(handles.Ia,'Value'); handles.Ia_txt.String = sprintf('Ia = %.2f',Ia);
    pn0= get(handles.Pn,'Value'); handles.Pn_txt.String = sprintf('Pn = %.2f',pn0);
end
SnnPD = pn0;
W0 = [1-pn0 pn0 0 0 0];
Tmax = 12;                            % change yearly
N =  600;

[SnPD, SnnPD, Ss, SrPD, SrnPD] = rk_Suicide_SD( W0);
Y = [SnPD;SnnPD;Ss;SrPD;SrnPD;Ss*c];
% TLTS = {'people with no suicidal thought \n but psychosocial distress',...
%         'people with no suicidal thought \n or psychosocial distress',...
%         'people with suicidal thought \n',...
%         'people in discovery program with \n psychosocial distress',...
%         'people in discovery program with \n no psychosocial distress',...
%         'people committed suicide \n'};
TLTS = {'S_n with PD',...
    'S_n without PD',...
    'S_s',...
    'S_r with PD',...
    'S_r without PD',...
    'S_{death}'};

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
    if k <5
        ylim([0,1])
        if k > 2
            ylim([0,0.4])
        end
    end
    if k ==5
        ylim([0,0.06])
    end
    xlim([0,12])
    
    
    xlabel('time (year)')
end


function [SnPD, SnnPD, Ss,SrPD, SrnPD] = rk_Suicide_SD(W)
%----------------------------------------------------------------------------------------------
% Runge-Kutta solver for the suicide ODEs
% W = Sn_PD Sn_nPD Ss Sr_Pd Sr_nPd
%----------------------------------------------------------------------------------------------
global Tmax N                                            % Tmax = maximum time (months), N = number of time points

n    = size(W,1);                                        % number of samples
Tmin = 0;                                                % minimum time
h    = (Tmax - Tmin) / N;                                % time increment
SnPD   = zeros(n,N);
SnnPD   = zeros(n,N); Ss = SnPD;
SrPD    = SnnPD;   SrnPD = SnnPD;                        % initialize the vecotors
for i = 1:N                                              % for each time step
    t       = Tmin + h*i;
    %     SnPD(:,i) = max(W(:,1), 0);                        % stock of non suicidal with PD
    %     SnnPD(:,i)= max(W(:,2), 0);                        % stock of suicidal without PD
    %     Ss(:,i)   = max(W(:,3), 0);                        % stock of suicidal
    %     SrPD(:,i) = max(W(:,4), 0);                        % stock of recovery with PD
    %     SrnPD(:,i)= max(W(:,5), 0);                        % stock of recovery without PD
    SnPD(:,i) = W(:,1);                        % stock of non suicidal with PD
    SnnPD(:,i)= W(:,2);                        % stock of suicidal without PD
    Ss(:,i)   = W(:,3);                        % stock of suicidal
    SrPD(:,i) = W(:,4);                        % stock of recovery with PD
    SrnPD(:,i)= W(:,5);                        % stock of recovery without PD
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
%   W     = Sn_PD Sn_nPD Ss Sr_Pd Sr_nPd
%----------------------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------------
global a0  ac ar aa b0 Ic Ia Ir c          % parameters for stocks into AT

SnPD    = W(:,1);                        % stock of non suicidal with PD
SnnPD   = W(:,2);                        % stock of suicidal without PD
Ss      = W(:,3);                        % stock of suicidal
SrPD    = W(:,4);                        % stock of recovery with PD
SrnPD   = W(:,5);                        % stock of recovery without PD


% rate constant
nnPD_nPD_rate      = a0;                 
nPD_nnPD_rate      = ac*Ic;
nPD_s_rate         = b0*(1 - aa*Ia);
s_rPD_rate         = a1*(1-ar*Ir);
s_rnPD_rate        = a1*((ar)*Ir);
rnPD_nnPD_rate     = 1;
rPD_nPD_rate       = 1;
s_suicide_rate     = c;

% flow
nnPD_nPD_flow      = nnPD_nPD_rate .* SnnPD;
nPD_nnPD_flow      = nPD_nnPD_rate .* SnPD ;
nPD_s_flow         = nPD_s_rate    .* SnPD;
s_rPD_flow         = s_rPD_rate    .* Ss;
s_rnPD_flow        = s_rnPD_rate   .* Ss;
rnPD_nnPD_flow     = rnPD_nnPD_rate.* SrnPD;
rPD_nPD_flow       = rPD_nPD_rate  .* SrPD;
s_suicide_flow     = s_suicide_rate.* Ss;
birth_flow         = s_suicide_flow;



% Derivatives
% temporary Derivatives
dSnnPD = -nnPD_nPD_flow + nPD_nnPD_flow + rnPD_nnPD_flow + birth_flow;            % dSnnPD/dt
dSnPD  = -nPD_nnPD_flow - nPD_s_flow + nnPD_nPD_flow +rPD_nPD_flow;
dSs    = -s_rPD_flow - s_rnPD_flow - s_suicide_flow + nPD_s_flow;
dSrPD  = - rPD_nPD_flow  + s_rPD_flow ;
dSrnPD = -rnPD_nnPD_flow + s_rnPD_flow;

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

% real Derivatives 
dSnnPD = -nnPD_nPD_flow + nPD_nnPD_flow + rnPD_nnPD_flow + birth_flow;            % dSnnPD/dt
dSnPD  = -nPD_nnPD_flow - nPD_s_flow + nnPD_nPD_flow +rPD_nPD_flow;
dSs    = -s_rPD_flow - s_rnPD_flow - s_suicide_flow + nPD_s_flow;
dSrPD  = - rPD_nPD_flow  + s_rPD_flow ;
dSrnPD = -rnPD_nnPD_flow + s_rnPD_flow;


dW(:,1) = dSnPD;
dW(:,2) = dSnnPD;
dW(:,3) = dSs;
dW(:,4) = dSrPD;
dW(:,5) = dSrnPD;
