%% GRAND_AVERAGE
% Grand average ERPs for Dual N-Back 
data         = 'C:\Users\sergi\OneDrive\Escritorio\ERP\CLEAN\BF';  
file     = '*.set';
chans = {'Fz','Cz','Pz'};
levels          = [0 1 2 3];
stimEventBase   = 'EVENT01';          
respEventType   = 'EVENT02';          
levelOrder = [0 1 2 3];       
epochWinSec     = [-0.5 2.5];          
baselineMs      = [-500 0];            
applySmoothing  = true;
smoothN         = 60;

% Low-pass filter
useLowpassERP   = true;    
lp.fc           = 3;       
lp.order        = 1;      
plotUnits       = 'µV';
lineWidth       = 1.8;


if saveFigs && ~exist(figOutDir, 'dir'); mkdir(figOutDir); end
d = dir(fullfile(data, file));
assert(~isempty(d), 'No .set files found in %s', data);

GA = struct();              
PerSub = struct();          
timeVecs = {};

if isempty(which('eeglab.m'))
    error('EEGLAB not found on path. addpath(genpath(''path/to/eeglab''))');
end

%% MAIN LOOP OVER SUBJECTS 
for s = 1:numel(d)
    setPath = fullfile(d(s).folder, d(s).name);
    fprintf('\n[%d/%d] Loading %s\n', s, numel(d), d(s).name);
    EEG = pop_loadset(setPath);

    for L = levels
        EEG_L = try_epoch_by_level(EEG, stimEventBase, L, epochWinSec);
        if isempty(EEG_L)
            EEG_L = try_epoch_by_block(EEG, respEventType, levelOrder, L, stimEventBase, epochWinSec);
        end
        if isempty(EEG_L)
            warning('No epochs for Level %d in %s. Skipping.', L, d(s).name);
            continue;
        end

        EEG_L = pop_rmbase(EEG_L, baselineMs);

        % Apply your 3 Hz low‑pass 
        if useLowpassERP
            EEG_L = applyLowpassERP(EEG_L, lp.fc, lp.order);
        end

        % ERP 
        erp = squeeze(mean(EEG_L.data, 3));
        if applySmoothing
            for c = 1:size(erp,1)
                erp(c,:) = smooth(erp(c,:), smoothN);
            end
        end

        % Store Fz/Cz/Pz 
        timeVec = linspace(epochWinSec(1)*1000, epochWinSec(2)*1000, size(erp,2)); % ms
        for ch = 1:numel(chans)
            chIdx = find(strcmpi({EEG_L.chanlocs.labels}, chans{ch}), 1);
            if isempty(chIdx), continue; end
            PerSub(s).level(L+1).(chans{ch}).time = timeVec;
            PerSub(s).level(L+1).(chans{ch}).erp  = erp(chIdx,:);
        end
        timeVecs{end+1} = timeVec; 
    end
end

%% ALIGN & GRAND AVERAGE 
assert(~isempty(timeVecs), 'No epochs gathered across subjects/levels.');

consTime = timeVecs{1};
for k = 2:numel(timeVecs)
    consTime = intersect(round(consTime,6), round(timeVecs{k},6));
end
assert(numel(consTime) >= 10, 'Too little overlap in time bases.');

for L = levels
    for ch = 1:numel(chans)
        chan = chans{ch};
        stack = [];
        for s = 1:numel(d)
            have = (isfield(PerSub, 'level') && numel(PerSub(s).level) >= (L+1) && ...
                    isfield(PerSub(s).level(L+1), chan) && ~isempty(PerSub(s).level(L+1).(chan)));
            if ~have, continue; end
            t = PerSub(s).level(L+1).(chan).time;
            y = PerSub(s).level(L+1).(chan).erp;
            y_al = interp1(t, y, consTime, 'linear', NaN);
            stack = [stack; y_al]; 
        end
        if isempty(stack), continue; end
        GA(L+1).(chan).time = consTime;
        GA(L+1).(chan).erp  = nanmean(stack, 1);
        GA(L+1).(chan).n    = sum(all(~isnan(stack),2));
    end
end

%% PLOTTING 
eventMarks = [0 500];   

for ch = 1:numel(chans)
    chan = chans{ch};
    f = figure('Color','w','Name',sprintf('GA - %s', chan));
    ax = axes('Parent', f); hold(ax, 'on');

    
    for L = levels
        if numel(GA) >= (L+1) && isfield(GA(L+1), chan) && ~isempty(GA(L+1).(chan))
            plot(ax, GA(L+1).(chan).time, GA(L+1).(chan).erp, 'LineWidth', lineWidth);
        else
            plot(ax, NaN,NaN, 'LineWidth', lineWidth); 
        end
    end

    xlabel(ax,'Time (ms)');
    ylabel(ax,plotUnits);

    h0 = yline(ax, 0, 'k-', 'LineWidth', 2.0);  hideFromLegend(h0);

    yl  = ylim(ax); yr = diff(yl);
    pad = 0.10 * yr; 
    ybot = yl(1) + pad;
    for em = eventMarks
        hx = xline(ax, em, 'k-', 'LineWidth', 2.0); hideFromLegend(hx);
        lbl = ternary(em==0, 'Audio', 'Visual');
        text(ax, em, ybot, lbl, 'Rotation', 90, 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold');
    end

    legend(ax, arrayfun(@(L) sprintf('Level %d', L), levels, 'UniformOutput', false), 'Location','Best');

    hold(ax, 'off');

    if saveFigs
        saveas(f, fullfile(figOutDir, sprintf('GA_%s_levels0-3.png', chan)));
    end
end


% Helpers

function EEG_L = try_epoch_by_level(EEG, stimBase, L, epochWinSec)
    EEG_L = [];
    evtypes = get_event_types(EEG); 
    target = [ ...
        sprintf('%s_L%d', stimBase, L) ...
        ; sprintf('%s_LEVEL%d', stimBase, L) ...
        ; string(sprintf('L%d', L)) ...
        ; string(sprintf('LEVEL%d', L)) ];
    if ~any(ismember(evtypes, target)), return; end
    firstType = target(find(ismember(target, evtypes), 1, 'first'));
    try
        EEG_L = pop_epoch(EEG, {char(firstType)}, epochWinSec);
    catch
        EEG_L = [];
    end
end

function EEG_L = try_epoch_by_block(EEG, respType, levelOrder, L, stimBase, epochWinSec)
    EEG_L = [];
    evtypes = get_event_types(EEG);
    hitResp = strcmpi(evtypes, string(respType));
    if ~any(hitResp), return; end
    respIdx = find(hitResp);
    if numel(respIdx) < 2, return; end
    respLat = arrayfun(@(k) EEG.event(k).latency, respIdx);

    blockStart = round(respLat(1:end-1));
    blockStop  = round(respLat(2:end) - 1);

    n = min([numel(levelOrder), numel(blockStart)]);
    if n < 1, return; end
    blockStart = blockStart(1:n);
    blockStop  = blockStop(1:n);
    levelOrder = levelOrder(1:n);

    idx = find(levelOrder == L, 1, 'first');
    if isempty(idx), return; end

    try
        EEG_block = pop_select(EEG, 'point', [blockStart(idx), blockStop(idx)]);
    catch
        EEG_block = [];
    end
    if isempty(EEG_block), return; end

    evtypesB = get_event_types(EEG_block);
    if ~any(strcmpi(evtypesB, string(stimBase)))
        return;
    end

    try
        EEG_L = pop_epoch(EEG_block, {stimBase}, epochWinSec);
    catch
        EEG_L = [];
    end
end

function evtypes = get_event_types(EEG)
    if ~isfield(EEG,'event') || isempty(EEG.event)
        evtypes = string.empty(1,0);
        return;
    end
    raw = {EEG.event.type};
    evtypes = strings(size(raw));
    for k = 1:numel(raw)
        v = raw{k};
        if isstring(v)
            evtypes(k) = v;
        elseif ischar(v)
            evtypes(k) = string(v);
        elseif isnumeric(v)
            evtypes(k) = string(num2str(v));
        else
            evtypes(k) = string("");
        end
    end
end

function Eout = applyLowpassERP(Ein, fc, orderBW)

    if isempty(Ein) || ~isfield(Ein,'data') || isempty(Ein.data)
        Eout = Ein; return;
    end
    fs = Ein.srate;
    [b,a] = butter(orderBW, fc/(fs/2), 'low');

    X = Ein.data;                       
    [~,~,N] = size(X);
    for n = 1:N
        xn = double(X(:,:,n));         
        X(:,:,n) = (filtfilt(b, a, xn')') ;  
    end
    Ein.data = X;
    Eout = Ein;
end

function hideFromLegend(h)

    try
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    catch
        
    end
end

function out = ternary(cond, a, b)

    if cond, out = a; else, out = b; end
end
