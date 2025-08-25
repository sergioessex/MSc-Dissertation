%% ERPs + TOPOPLOTS
% ERPs at Level 0-3 (Fz,Cz,Pz)
set(groot,'defaultAxesFontName','Arial','defaultAxesFontSize',10, ...
          'defaultAxesLineWidth',0.9,'defaultLineLineWidth',2.0);
set(groot,'defaultFigureColor','w');

% Parameters 
epochWinSec = [-0.5 2.5];      
baselineMs  = [-500 0];        
eventMarks  = [0 500];        
chansERP    = {'Fz','Cz','Pz'};
smoothN     = 60;              
p3SearchMs  = [800 1000];  p3HalfWin = 10;  p3RefChan = 'Pz';
n2SearchMs  = [300 500];   n2HalfWin = 10;  n2RefChan = 'Cz';
erpYLims    = [-4 4];
levelNames  = {'Level 1','Level 2','Level 3'};
levelColors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880];

% Low pass filter 
lp.fc    = 3;       
lp.order = 1;      

% Helper
getChan = @(E,lab) find(strcmpi({E.chanlocs.labels}, lab), 1);

% Time ranges per level (BLOCK LATENCY)
Block_latency = [EEG.event(ismember({EEG.event.type},'EVENT02')).latency];
last_event    = EEG.event(end).latency + 3*500;
ranges_pts = { ...
    [Block_latency(1),  Block_latency(2)-1];  
    [Block_latency(2),  Block_latency(6)-1];   
    [Block_latency(6),  Block_latency(10)-1];  
    [Block_latency(10), last_event]};         

% Build level datasets
EEGs = cell(1,4);
for L = 1:4
    E = pop_select(EEG,'point',ranges_pts{L});
    E = pop_epoch(E, {'EVENT01'}, epochWinSec);
    E = pop_rmbase(E, baselineMs);

    % 3 Hz low-pass for ERPs 
    E = applyLowpassERP(E, lp.fc, lp.order);

    EEGs{L} = E;
end
EEG_0 = EEGs{1}; EEG_1 = EEGs{2}; EEG_2 = EEGs{3}; EEG_3 = EEGs{4};
levels = {EEG_1, EEG_2, EEG_3};
t = EEG_1.times(:)';           

% ERPs per channel and level 
ERP = cell(numel(chansERP),3);
for c = 1:numel(chansERP)
    for L = 1:3
        E = levels{L};
        ch = getChan(E, chansERP{c});
        assert(~isempty(ch), 'Channel %s not found in dataset.', chansERP{c});
        data = squeeze(E.data(ch,:,:));
        erp  = mean(data,2);
        ERP{c,L} = normalize(smooth(erp, smoothN));
    end
end

% ERPs Level 0 per channel (Fz, Cz, Pz) 
ERP0 = cell(1,numel(chansERP));
for c = 1:numel(chansERP)
    ch0 = getChan(EEG_0, chansERP{c});
    assert(~isempty(ch0), 'Channel %s not found in Level 0.', chansERP{c});
    data0 = squeeze(EEG_0.data(ch0,:,:));      
    erp0  = mean(data0,2);
    ERP0{c} = normalize(smooth(erp0, smoothN));
end

% Peaks at P300 and N200
maskP3 = (t >= p3SearchMs(1)) & (t <= p3SearchMs(2));
maskN2 = (t >= n2SearchMs(1)) & (t <= n2SearchMs(2));
pzRow  = find(strcmpi(chansERP,p3RefChan),1);
czRow  = find(strcmpi(chansERP,n2RefChan),1);

tPkP3 = zeros(1,3); tPkN2 = zeros(1,3);
for L = 1:3
    [~, iMax] = max(ERP{pzRow,L}(maskP3));  tt = t(maskP3);  tPkP3(L) = tt(iMax);
    [~, iMin] = min(ERP{czRow,L}(maskN2));  tt = t(maskN2);  tPkN2(L) = tt(iMin);
end

% Peaks for Level 0 
[~, iMax0] = max(ERP0{pzRow}(maskP3));  ttP3 = t(maskP3);  tPkP3_0 = ttP3(iMax0);
[~, iMin0] = min(ERP0{czRow}(maskN2));  ttN2 = t(maskN2);  tPkN2_0 = ttN2(iMin0);

% Topomap vectors (+-10 ms around peaks)
TopoP3 = cell(1,3);  TopoN2 = cell(1,3);
for L = 1:3
    E = levels{L}; erpAll = squeeze(mean(E.data,3)); 
    mP3 = (t >= tPkP3(L)-p3HalfWin) & (t <= tPkP3(L)+p3HalfWin);
    mN2 = (t >= tPkN2(L)-n2HalfWin) & (t <= tPkN2(L)+n2HalfWin);
    TopoP3{L} = mean(erpAll(:,mP3),2);
    TopoN2{L} = mean(erpAll(:,mN2),2);
end

erpAll0 = squeeze(mean(EEG_0.data,3));
mP3_0 = (t >= tPkP3_0-p3HalfWin) & (t <= tPkP3_0+p3HalfWin);
mN2_0 = (t >= tPkN2_0-n2HalfWin) & (t <= tPkN2_0+n2HalfWin);
TopoP3_0 = mean(erpAll0(:,mP3_0),2);
TopoN2_0 = mean(erpAll0(:,mN2_0),2);
allP3 = [TopoP3_0; TopoP3{1}; TopoP3{2}; TopoP3{3}];
allN2 = [TopoN2_0; TopoN2{1}; TopoN2{2}; TopoN2{3}];
climsP3 = [-max(abs(allP3)),  max(abs(allP3))];
climsN2 = [-max(abs(allN2)),  max(abs(allN2))];


% ERP figures: one per channel (Fz, Cz, Pz) 

ERPfigure('Fz','ERPs_Fz');
ERPfigure('Cz','ERPs_Cz');
ERPfigure('Pz','ERPs_Pz');
makeERPfigure_Level0_channels('ERP_Level0','ERPs_Level0');

% Topoplots
figTopo = figure('Name','Topoplots','Color','w','Renderer','painters');
tl2 = tiledlayout(2,4,'TileSpacing','compact','Padding','compact'); 

smallFS = 9;   
labelY  = 1.20;

% P300 row 
% Level 0
ax = nexttile; set(ax,'Color','w');
topoplot(TopoP3_0, EEG_0.chanlocs, 'maplimits', climsP3, ...
         'electrodes','on','style','map','shading','interp');
text(ax, 0.5, labelY, 'LEVEL 0 — P300', ...
     'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top', ...
     'FontSize', smallFS, 'FontWeight', 'bold','Clipping','on');

for L = 1:3
    ax = nexttile; set(ax,'Color','w');
    topoplot(TopoP3{L}, levels{L}.chanlocs, 'maplimits', climsP3, ...
             'electrodes','on','style','map','shading','interp');
    text(ax, 0.5, labelY, sprintf('LEVEL %d — P300', L), ...
         'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top', ...
         'FontSize', smallFS, 'FontWeight', 'bold','Clipping','on');
end
colormap(figTopo, parula);
cbP3 = colorbar; cbP3.Location = 'eastoutside';
cbP3.Label.String = sprintf('Mean \\muV (P300, \\pm%d ms)', p3HalfWin);

% N200 row 
% Level 0
ax = nexttile; set(ax,'Color','w');
topoplot(TopoN2_0, EEG_0.chanlocs, 'maplimits', climsN2, ...
         'electrodes','on','style','map','shading','interp');
text(ax, 0.5, labelY, 'LEVEL 0 — N200', ...
     'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top', ...
     'FontSize', smallFS, 'FontWeight', 'bold','Clipping','on');

for L = 1:3
    ax = nexttile; set(ax,'Color','w');
    topoplot(TopoN2{L}, levels{L}.chanlocs, 'maplimits', climsN2, ...
             'electrodes','on','style','map','shading','interp');
    text(ax, 0.5, labelY, sprintf('LEVEL %d — N200', L), ...
         'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top', ...
         'FontSize', smallFS, 'FontWeight', 'bold','Clipping','on');
end
cbN2 = colorbar; cbN2.Location = 'eastoutside';
cbN2.Label.String = sprintf('Mean \\muV (N200, \\pm%d ms)', n2HalfWin);

%% Local functions 
function E = applyLowpassERP(E, fc, order)

    if isempty(E) || ~isfield(E,'data') || isempty(E.data)
        return; 
    end

    fs = E.srate;
    [b,a] = butter(order, fc/(fs/2), 'low');

    xOrigClass = class(E.data);
    X = double(E.data);                
    [C,T,N] = size(X);

    % Vectorize across channels*trials.
    X2  = permute(X, [2 1 3]);         
    X2  = reshape(X2, T, C*N);         
    X2f = filtfilt(b, a, X2);          
    Xf  = reshape(X2f, T, C, N);       
    Xf  = permute(Xf, [2 1 3]);        

    E.data = cast(Xf, xOrigClass);
end

function ERPfigure(chanLabel, baseFileName)

    t            = evalin('base','t');
    ERP          = evalin('base','ERP');
    chansERP     = evalin('base','chansERP');
    erpYLims     = evalin('base','erpYLims');
    levelColors  = evalin('base','levelColors');
    levelNames   = evalin('base','levelNames');
    n2SearchMs   = evalin('base','n2SearchMs');
    p3SearchMs   = evalin('base','p3SearchMs');
    eventMarks   = evalin('base','eventMarks');
    epochWinSec  = evalin('base','epochWinSec');
    tPkP3        = evalin('base','tPkP3');
    tPkN2        = evalin('base','tPkN2');

    c = find(strcmpi(chansERP,chanLabel),1);

    fig = figure('Name',['ERP_' chanLabel],'Color','w','Renderer','painters');
    tl  = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
    ax  = nexttile; hold(ax,'on'); set(ax,'Color','w');

    yl  = erpYLims; yr = diff(yl);

    pN2 = patch(ax,[n2SearchMs(1) n2SearchMs(2) n2SearchMs(2) n2SearchMs(1)], ...
                   [yl(1) yl(1) yl(2) yl(2)],[0.80 0.90 1.00], ...
                   'EdgeColor','none','FaceAlpha',0.25);
    pP3 = patch(ax,[p3SearchMs(1) p3SearchMs(2) p3SearchMs(2) p3SearchMs(1)], ...
                   [yl(1) yl(1) yl(2) yl(2)],[1.00, 0.80, 0.80], ...
                   'EdgeColor','none','FaceAlpha',0.25);
    hideFromLegend(pN2); hideFromLegend(pP3);

    for L = 1:3
        plot(ax, t, ERP{c,L}, 'Color', levelColors(L,:), ...
            'DisplayName', levelNames{L}, 'LineWidth', 2.2);
    end

    h0 = yline(ax,0,'k-','LineWidth',2.0);  hideFromLegend(h0);
    pad = 0.10*yr; ybot = yl(1)+pad;
    for em = eventMarks
        hx = xline(ax,em,'k-','LineWidth',2.0); hideFromLegend(hx);
        lbl = ternary(em==0,'Audio','Visual');
        text(ax,em,ybot,lbl,'Rotation',90,'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom','FontSize',10,'FontWeight','bold');
    end

    for L = 1:3
        [~,iN2] = min(abs(t - tPkN2(L))); yN2 = ERP{c,L}(iN2);
        plot(ax,tPkN2(L),yN2,'v','MarkerSize',7,'MarkerFaceColor',levelColors(L,:),'MarkerEdgeColor','k');
        [~,iP3] = min(abs(t - tPkP3(L))); yP3 = ERP{c,L}(iP3);
        plot(ax,tPkP3(L),yP3,'^','MarkerSize',7,'MarkerFaceColor',levelColors(L,:),'MarkerEdgeColor','k');
    end

    ylbl = yl(2) - 0.08*yr;
    text(ax,mean(n2SearchMs), ylbl, 'N2', 'HorizontalAlignment','center', ...
         'VerticalAlignment','top','FontSize',11,'FontWeight','bold');
    text(ax,mean(p3SearchMs), ylbl, 'P3', 'HorizontalAlignment','center', ...
         'VerticalAlignment','top','FontSize',11,'FontWeight','bold');

    grid(ax,'off'); box(ax,'off');
    xlim(ax,[epochWinSec*1000]); ylim(ax,yl);
    ylabel(ax,'\muV'); xlabel(ax,'Time [ms]');

    lgd = legend(ax, {'Level 1','Level 2','Level 3'}, ...
        'Location','northeast','Orientation','vertical','Box','on', ...
        'FontSize',12,'FontWeight','bold','TextColor','k');
    lgd.Color     = 'w';
    lgd.EdgeColor = [0 0 0];
    lgd.LineWidth = 0.75;
    try, lgd.ItemTokenSize = [40,18]; end
end

function makeERPfigure_Level0_channels(figName, baseFileName)
    t           = evalin('base','t');
    ERP0        = evalin('base','ERP0');           
    chansERP    = evalin('base','chansERP');
    erpYLims    = evalin('base','erpYLims');
    levelColors = evalin('base','levelColors');    
    n2SearchMs  = evalin('base','n2SearchMs');
    p3SearchMs  = evalin('base','p3SearchMs');
    eventMarks  = evalin('base','eventMarks');
    epochWinSec = evalin('base','epochWinSec');
    tPkP3_0     = evalin('base','tPkP3_0');
    tPkN2_0     = evalin('base','tPkN2_0');

    fig = figure('Name',figName,'Color','w','Renderer','painters');
    tl  = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
    ax  = nexttile; hold(ax,'on'); set(ax,'Color','w');

    yl  = erpYLims; yr = diff(yl);

    pN2 = patch(ax,[n2SearchMs(1) n2SearchMs(2) n2SearchMs(2) n2SearchMs(1)], ...
                   [yl(1) yl(1) yl(2) yl(2)],[0.80 0.90 1.00], ...
                   'EdgeColor','none','FaceAlpha',0.25);
    pP3 = patch(ax,[p3SearchMs(1) p3SearchMs(2) p3SearchMs(2) p3SearchMs(1)], ...
                   [yl(1) yl(1) yl(2) yl(2)],[1.00, 0.80, 0.80], ...
                   'EdgeColor','none','FaceAlpha',0.25);
    hideFromLegend(pN2); hideFromLegend(pP3);

    for c = 1:numel(chansERP)
        plot(ax, t, ERP0{c}, 'Color', levelColors(c,:), ...
            'DisplayName', chansERP{c}, 'LineWidth', 2.2);
    end

    h0 = yline(ax,0,'k-','LineWidth',2.0);  hideFromLegend(h0);
    pad = 0.10*yr; ybot = yl(1)+pad;
    for em = eventMarks
        hx = xline(ax,em,'k-','LineWidth',2.0); hideFromLegend(hx);
        lbl = ternary(em==0,'Audio','Visual');
        text(ax,em,ybot,lbl,'Rotation',90,'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom','FontSize',10,'FontWeight','bold');
    end

    for c = 1:numel(chansERP)
        [~,iN2c] = min(abs(t - tPkN2_0)); yN2c = ERP0{c}(iN2c);
        plot(ax,tPkN2_0,yN2c,'v','MarkerSize',7,'MarkerFaceColor',levelColors(c,:),'MarkerEdgeColor','k');
        [~,iP3c] = min(abs(t - tPkP3_0)); yP3c = ERP0{c}(iP3c);
        plot(ax,tPkP3_0,yP3c,'^','MarkerSize',7,'MarkerFaceColor',levelColors(c,:),'MarkerEdgeColor','k');
    end

    ylbl = yl(2) - 0.08*yr;
    text(ax,mean(n2SearchMs), ylbl, 'N2', 'HorizontalAlignment','center', ...
         'VerticalAlignment','top','FontSize',11,'FontWeight','bold');
    text(ax,mean(p3SearchMs), ylbl, 'P3', 'HorizontalAlignment','center', ...
         'VerticalAlignment','top','FontSize',11,'FontWeight','bold');

    grid(ax,'off'); box(ax,'off');
    xlim(ax,[epochWinSec*1000]); ylim(ax,yl);
    ylabel(ax,'\muV'); xlabel(ax,'Time [ms]');

    lgd = legend(ax, chansERP, ...
        'Location','northeast','Orientation','vertical','Box','on', ...
        'FontSize',12,'FontWeight','bold','TextColor','k');
    lgd.Color     = 'w';
    lgd.EdgeColor = [0 0 0];
    lgd.LineWidth = 0.75;
    try, lgd.ItemTokenSize = [40,18]; end
end

function hideFromLegend(h)
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

function out = ternary(cond,a,b)
    if cond, out=a; else, out=b; end
end

