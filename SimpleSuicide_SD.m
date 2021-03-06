function SimpleSuicide_SD(handles)
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
a0 = 1;
ar = 0.5;
ac = 0.5;
aa = 0.5;
b0 = 0.2;

Ir = 0.4;
Ic = 0.0;
Ia = 0.0;

pn0= 0.3;

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

Is0   = {Ir,Ic,Ia};
Ilbl  = {'recovery program','community','awareness','community'};
nplot = 30;
Iplotx= linspace(0,1,nplot);
Iploty= linspace(0,1,nplot+1);

if ~nargin
    fig = figure(100); fig.Name = 'ES'; fig.Color = 'w';
    for k=1:3
        Is     = Is0;
        Is{k}  = Iplotx;
        [ pntots{k},Sn,EStots{k},ESr,ESc,ESa,pntot0s{k},EStot0s{k},ESr0,ESc0,ESa0,ES0 ] = fES(a0,ar,ac,aa,b0,Is,pn0);
        
        subplot(3,2,k*2-1)
        plot(Iplotx,(EStot0s{k})/ES0*100,'r-',Iplotx,(EStots{k})/ES0*100,'b-','LineWidth',2)        
        ax = gca; ax.FontSize = 14; ax.XLim = [0 1];
        xlabel([Ilbl{k} ' intervention effort'])
        ylabel('% reduction')
        if k==1
            title('Model comparison: effect size on suicide rate')
            legend('linear model','SD model','Location','SE')
        end
        
        subplot(3,2,k*2)
        plot(Iplotx,(pn0-pntot0s{k})/pn0*100,'r-',Iplotx,(pn0-pntots{k})/pn0*100,'b-','LineWidth',2)
        ax = gca; ax.FontSize = 14; ax.XLim = [0 1];
        xlabel([Ilbl{k} ' intervention effort'])
        ylabel('% reduction')
        if k==1
            title('Model comparison: effect size on PD')
            legend('linear model','SD model','Location','SE')
        end
    end
    
    k = 4; % do two changes in interventions at the same time
    Is = Is0; Is{2} = Iplotx; Is{3} = Iplotx(end:-1:1);
    [ pntots{k},Sn,EStots{k},ESr,ESc,ESa,pntot0s{k},EStot0s{k},ESr0,ESc0,ESa0,ES0 ] = fES(a0,ar,ac,aa,b0,Is,pn0);

    fig = figure(101); fig.Name = 'ES two Is'; fig.Color = 'w';
    plot(Iplotx,(EStot0s{k})/ES0*100,'r-',Iplotx,(EStots{k})/ES0*100,'b-','LineWidth',2)
    ax = gca; ax.FontSize = 14; ax.XLim = [0 1]; ax.YLim = [0 57];
    xlabel([Ilbl{k} ' intervention effort'])
    ylabel('effect size on suicide rate (% reduction)')
    
    fig = figure(102); fig.Name = 'Pn two Is'; fig.Color = 'w';
    plot(Iplotx,(pn0-pntot0s{k})/pn0*100,'r-',Iplotx,(pn0-pntots{k})/pn0*100,'b-','LineWidth',2)
    ax = gca; ax.FontSize = 14; ax.XLim = [0 1]; ax.YLim = [0 50]; ax.YTick = 0:10:50;
    ylabel('effect size on PD (% reduction)')
    xlabel([Ilbl{k} ' intervention effort'])

    keyboard
    ax.XDir = 'reverse'; xlabel('awareness intervention effort')
    
    return
end

[ Ixy{1},Ixy{2} ] = meshgrid(Iplotx,Iploty);
if nargin
    h_axes = {handles.axes1,handles.axes2,handles.axes3,...
        handles.axes4,handles.axes5,handles.axes6,handles.axes7,handles.axes8,handles.axes9};
else
    fig = figure(110); fig.Name = 'ES 3D'; clf, h_axes = [];
end

for k=1:3
    Is   = Is0;
    notk = setdiff(1:3,k);
    for j=1:2
        Is{notk(j)} = Ixy{j};
    end
    [ pntot,Sn,EStot,ESr,ESc,ESa,pntot0,EStot0,ESr0,ESc0,ESa0,ES0 ] = fES(a0,ar,ac,aa,b0,Is,pn0);
    
    if k<3, zmax = max(EStot(:)); end
    plotImagesc(k,1,Iplotx,Iploty,EStot, Ilbl{notk(1)},Ilbl{notk(2)},Ilbl{k},Is{k},h_axes,nargin,zmax)
    plotImagesc(k,2,Iplotx,Iploty,EStot0,Ilbl{notk(1)},Ilbl{notk(2)},Ilbl{k},Is{k},h_axes,nargin,zmax)
    plotImagesc(k,3,Iplotx,Iploty,pntot, Ilbl{notk(1)},Ilbl{notk(2)},Ilbl{k},Is{k},h_axes,nargin,pn0 )
end


function plotImagesc(krow,kcol,Iplotx,Iploty,Z,Ilblx,Ilbly,Ilblk,Isk,h_axes,Nargin,zmax)
%----------------------------------------------------------------------------------------------
% plot imagesc of Z
%----------------------------------------------------------------------------------------------
if Nargin, axes(h_axes{krow+kcol*3-3}), else, subplot(3,3,sub2ind([3 3],kcol,krow)), end
imagesc(Iplotx,Iploty,Z)
if ~isempty(zmax), caxis([0 zmax]), end
xlabel(Ilblx), ylabel(Ilbly), title(sprintf('%s, I=%.2f',Ilblk,Isk))
ax = gca; ax.YDir = 'normal';


function [ pntot,Sn,EStot,ESr,ESc,ESa,pntot0,EStot0,ESr0,ESc0,ESa0,ES0 ] = fES(a0,ar,ac,aa,b0,Is,pn0)
%----------------------------------------------------------------------------------------------
% function for ES and pn and Sn
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
Ir   = Is{1}; Ic = Is{2}; Ia = Is{3};
b    = @(ia) b0*(1 - aa*ia);
pn   = @(ia,ir,ic) a0*pn0 ./ (a0 + ac*ic + ar*ir.*b(ia));
pntot= pn(Ia,Ir,Ic);
btot = b(Ia);
Sn   = 1./(1 + 2*btot.*pntot);

pntot0 = (1-ac*Ic)*pn0; % pn due to community intervention
ES0   = b0*        pn0;  % baseline rate coefficient for flow into the suicidal stock
ESa0  = b0*  aa*Ia*pn0;  % decrease in the coefficient due to awareness
ESc0  = b0*  ac*Ic*pn0;  %                                    community
ESr0  = b0^2*ar*Ir*pn0;  %                                    recovery program
EStot0 = ESr0 + ESc0 + ESa0;

ES   = @(ia,ic,ir) (ES0 - b(ia).*pn(ia,ir,ic));

EStot= ES(Ia,Ic,Ir);
ESr  = ES(0, 0, Ir);
ESa  = ES(Ia,0, 0 );
ESc  = ES(0, Ic,0 );