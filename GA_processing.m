dataDir = 'C:\Users\sergi\OneDrive\Escritorio\ERP\CLEAN\BF';
files = dir(fullfile(dataDir,'*.set'));

N = numel(files);
t = [];
GA_Fz0 = []; GA_Fz1 = []; GA_Fz2 = []; GA_Fz3 = [];
GA_Cz0 = []; GA_Cz1 = []; GA_Cz2 = []; GA_Cz3 = [];
GA_Pz0 = []; GA_Pz1 = []; GA_Pz2 = []; GA_Pz3 = [];

for s = 1:N
    EEG = pop_loadset(fullfile(files(s).folder, files(s).name));
    Block_latency = [EEG.event(ismember({EEG.event.type},'EVENT02')).latency];
    last_event = EEG.event(end).latency + 3*EEG.srate;

    EEG_0 = pop_select(EEG,'point',[Block_latency(1),Block_latency(2)-1]);
    EEG_0 = pop_epoch(EEG_0,{'EVENT01'},[-0.5 2.5]);
    EEG_0 = pop_rmbase(EEG_0,[-500 0]);
    ERP_0 = squeeze(mean(EEG_0.data,3));
    ERP_0_Fz = smooth(ERP_0(ismember({EEG_0.chanlocs.labels},'Fz'),:),70);
    ERP_0_Cz = smooth(ERP_0(ismember({EEG_0.chanlocs.labels},'Cz'),:),70);
    ERP_0_Pz = smooth(ERP_0(ismember({EEG_0.chanlocs.labels},'Pz'),:),70);

    EEG_1 = pop_select(EEG,'point',[Block_latency(2),Block_latency(6)-1]);
    EEG_1 = pop_epoch(EEG_1,{'EVENT01'},[-0.5 2.5]);
    EEG_1 = pop_rmbase(EEG_1,[-500 0]);
    ERP_1 = squeeze(mean(EEG_1.data,3));
    ERP_1_Fz = smooth(ERP_1(ismember({EEG_1.chanlocs.labels},'Fz'),:),70);
    ERP_1_Cz = smooth(ERP_1(ismember({EEG_1.chanlocs.labels},'Cz'),:),70);
    ERP_1_Pz = smooth(ERP_1(ismember({EEG_1.chanlocs.labels},'Pz'),:),70);

    EEG_2 = pop_select(EEG,'point',[Block_latency(6),Block_latency(10)-1]);
    EEG_2 = pop_epoch(EEG_2,{'EVENT01'},[-0.5 2.5]);
    EEG_2 = pop_rmbase(EEG_2,[-500 0]);
    ERP_2 = squeeze(mean(EEG_2.data,3));
    ERP_2_Fz = smooth(ERP_2(ismember({EEG_2.chanlocs.labels},'Fz'),:),70);
    ERP_2_Cz = smooth(ERP_2(ismember({EEG_2.chanlocs.labels},'Cz'),:),70);
    ERP_2_Pz = smooth(ERP_2(ismember({EEG_2.chanlocs.labels},'Pz'),:),70);

    EEG_3 = pop_select(EEG,'point',[Block_latency(10),last_event]);
    EEG_3 = pop_epoch(EEG_3,{'EVENT01'},[-0.5 2.5]);
    EEG_3 = pop_rmbase(EEG_3,[-500 0]);
    ERP_3 = squeeze(mean(EEG_3.data,3));
    ERP_3_Fz = smooth(ERP_3(ismember({EEG_3.chanlocs.labels},'Fz'),:),70);
    ERP_3_Cz = smooth(ERP_3(ismember({EEG_3.chanlocs.labels},'Cz'),:),70);
    ERP_3_Pz = smooth(ERP_3(ismember({EEG_3.chanlocs.labels},'Pz'),:),70);

    if isempty(t); t = EEG_0.times; end

    if s==1
        GA_Fz0 = ERP_0_Fz; GA_Fz1 = ERP_1_Fz; GA_Fz2 = ERP_2_Fz; GA_Fz3 = ERP_3_Fz;
        GA_Cz0 = ERP_0_Cz; GA_Cz1 = ERP_1_Cz; GA_Cz2 = ERP_2_Cz; GA_Cz3 = ERP_3_Cz;
        GA_Pz0 = ERP_0_Pz; GA_Pz1 = ERP_1_Pz; GA_Pz2 = ERP_2_Pz; GA_Pz3 = ERP_3_Pz;
    else
        GA_Fz0 = GA_Fz0 + ERP_0_Fz; GA_Fz1 = GA_Fz1 + ERP_1_Fz; GA_Fz2 = GA_Fz2 + ERP_2_Fz; GA_Fz3 = GA_Fz3 + ERP_3_Fz;
        GA_Cz0 = GA_Cz0 + ERP_0_Cz; GA_Cz1 = GA_Cz1 + ERP_1_Cz; GA_Cz2 = GA_Cz2 + ERP_2_Cz; GA_Cz3 = GA_Cz3 + ERP_3_Cz;
        GA_Pz0 = GA_Pz0 + ERP_0_Pz; GA_Pz1 = GA_Pz1 + ERP_1_Pz; GA_Pz2 = GA_Pz2 + ERP_2_Pz; GA_Pz3 = GA_Pz3 + ERP_3_Pz;
    end
end

GA_Fz0 = GA_Fz0 / N; GA_Fz1 = GA_Fz1 / N; GA_Fz2 = GA_Fz2 / N; GA_Fz3 = GA_Fz3 / N;
GA_Cz0 = GA_Cz0 / N; GA_Cz1 = GA_Cz1 / N; GA_Cz2 = GA_Cz2 / N; GA_Cz3 = GA_Cz3 / N;
GA_Pz0 = GA_Pz0 / N; GA_Pz1 = GA_Pz1 / N; GA_Pz2 = GA_Pz2 / N; GA_Pz3 = GA_Pz3 / N;

figure; hold on
plot(t, GA_Fz0); plot(t, GA_Fz1); plot(t, GA_Fz2); plot(t, GA_Fz3);
xline(0,'k'); yline(0,'k'); xline(500,'k'); legend('Level 0','Level 1','Level 2','Level 3'); title('Fz'); xlabel('Time (ms)'); ylabel('\muV')

figure; hold on
plot(t, GA_Cz0); plot(t, GA_Cz1); plot(t, GA_Cz2); plot(t, GA_Cz3);
xline(0,'k'); yline(0,'k'); xline(500,'k'); legend('Level 0','Level 1','Level 2','Level 3'); title('Cz'); xlabel('Time (ms)'); ylabel('\muV')

figure; hold on
plot(t, GA_Pz0); plot(t, GA_Pz1); plot(t, GA_Pz2); plot(t, GA_Pz3);
xline(0,'k'); yline(0,'k'); xline(500,'k'); legend('Level 0','Level 1','Level 2','Level 3'); title('Pz'); xlabel('Time (ms)'); ylabel('\muV')
