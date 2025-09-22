Block_latency = [EEG.event(ismember({EEG.event.type},'EVENT02')).latency];
last_event = EEG.event(end).latency + 3*EEG.srate;

EEG_0 = pop_select(EEG,'point',[Block_latency(1),Block_latency(2)-1]);
EEG_0 = pop_epoch(EEG_0,{'EVENT01'},[-0.5 2.5]);
EEG_0 = pop_rmbase(EEG_0,[-500 0]);
ERP_0 = squeeze(mean(EEG_0.data,3));
ERP_0_Fz = smooth(ERP_0(ismember({EEG_0.chanlocs.labels},'Fz'),:),70);
ERP_0_Cz = smooth(ERP_0(ismember({EEG_0.chanlocs.labels},'Cz'),:),70);
ERP_0_Pz = smooth(ERP_0(ismember({EEG_0.chanlocs.labels},'Pz'),:),70);
figure; hold on
plot(EEG_0.times,ERP_0_Fz)
plot(EEG_0.times,ERP_0_Cz)
plot(EEG_0.times,ERP_0_Pz)
xline(0,'k'); yline(0,'k'); xline(500,'k')
legend('Fz','Cz','Pz'); title('Level 0')

EEG_1 = pop_select(EEG,'point',[Block_latency(2),Block_latency(6)-1]);
EEG_1 = pop_epoch(EEG_1,{'EVENT01'},[-0.5 2.5]);
EEG_1 = pop_rmbase(EEG_1,[-500 0]);
ERP_1 = squeeze(mean(EEG_1.data,3));
ERP_1_Fz = smooth(ERP_1(ismember({EEG_1.chanlocs.labels},'Fz'),:),70);
ERP_1_Cz = smooth(ERP_1(ismember({EEG_1.chanlocs.labels},'Cz'),:),70);
ERP_1_Pz = smooth(ERP_1(ismember({EEG_1.chanlocs.labels},'Pz'),:),70);
figure; hold on
plot(EEG_1.times,ERP_1_Fz)
plot(EEG_1.times,ERP_1_Cz)
plot(EEG_1.times,ERP_1_Pz)
xline(0,'k'); yline(0,'k'); xline(500,'k')
legend('Fz','Cz','Pz'); title('Level 1')

EEG_2 = pop_select(EEG,'point',[Block_latency(6),Block_latency(10)-1]);
EEG_2 = pop_epoch(EEG_2,{'EVENT01'},[-0.5 2.5]);
EEG_2 = pop_rmbase(EEG_2,[-500 0]);
ERP_2 = squeeze(mean(EEG_2.data,3));
ERP_2_Fz = smooth(ERP_2(ismember({EEG_2.chanlocs.labels},'Fz'),:),70);
ERP_2_Cz = smooth(ERP_2(ismember({EEG_2.chanlocs.labels},'Cz'),:),70);
ERP_2_Pz = smooth(ERP_2(ismember({EEG_2.chanlocs.labels},'Pz'),:),70);
figure; hold on
plot(EEG_2.times,ERP_2_Fz)
plot(EEG_2.times,ERP_2_Cz)
plot(EEG_2.times,ERP_2_Pz)
xline(0,'k'); yline(0,'k'); xline(500,'k')
legend('Fz','Cz','Pz'); title('Level 2')

EEG_3 = pop_select(EEG,'point',[Block_latency(10),last_event]);
EEG_3 = pop_epoch(EEG_3,{'EVENT01'},[-0.5 2.5]);
EEG_3 = pop_rmbase(EEG_3,[-500 0]);
ERP_3 = squeeze(mean(EEG_3.data,3));
ERP_3_Fz = smooth(ERP_3(ismember({EEG_3.chanlocs.labels},'Fz'),:),70);
ERP_3_Cz = smooth(ERP_3(ismember({EEG_3.chanlocs.labels},'Cz'),:),70);
ERP_3_Pz = smooth(ERP_3(ismember({EEG_3.chanlocs.labels},'Pz'),:),70);
figure; hold on
plot(EEG_3.times,ERP_3_Fz)
plot(EEG_3.times,ERP_3_Cz)
plot(EEG_3.times,ERP_3_Pz)
xline(0,'k'); yline(0,'k'); xline(500,'k')
legend('Fz','Cz','Pz'); title('Level 3')
