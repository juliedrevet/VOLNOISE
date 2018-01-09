function verification_jitters(expe)

% % % 
% % % iti   = [.5 2.2 3.9]; % inter-trial interval (fixation point before stim onset)
% % % stout = [.4 2.1 3.8]-.3;

l1 = 1;
l2 = 1;
l3 = 1;
for iblck = 1:3
    for ivol = 1:3
        idx = find(expe.blck(iblck).switch_seq(((ivol-1)*60+1):ivol*60));
        idx = idx(idx>3);
        idx = idx(idx<59);
% % %         itis = expe.clck(iblck).tstim-expe.clck(iblck).tfix1(1:end-1);
% % %         stouts = expe.clck(iblck).toutc-expe.clck(iblck).tstim-1;
% % %         iti_jitts   = itis(((ivol-1)*60+1):ivol*60);
% % %         stout_jitts = stouts(((ivol-1)*60+1):ivol*60);
        iti_jitts   = expe.stim(iblck).iti_jitters(((ivol-1)*60+1):ivol*60);
        stout_jitts = expe.stim(iblck).stout_jitters(((ivol-1)*60+1):ivol*60);
        if expe.blck(iblck).volperm(ivol)==1
            
            for i = 1:length(idx)
                iti.vol1(l1,:) = iti_jitts((idx(i)-3):(idx(i)+2));
                stout.vol1(l1,:) = stout_jitts((idx(i)-3):(idx(i)+2));
                l1 = l1+1;
            end
        elseif expe.blck(iblck).volperm(ivol)==2
            for i = 1:length(idx)
                iti.vol2(l2,:) = iti_jitts((idx(i)-3):(idx(i)+2));
                stout.vol2(l2,:) = stout_jitts((idx(i)-3):(idx(i)+2));
                l2 = l2+1;
            end
        elseif expe.blck(iblck).volperm(ivol)==3
            for i = 1:length(idx)
                iti.vol3(l3,:) = iti_jitts((idx(i)-3):(idx(i)+2));
                stout.vol3(l3,:) = stout_jitts((idx(i)-3):(idx(i)+2));
                l3 = l3+1;
            end
        end
    end
end

iti.vol1(iti.vol1<1) = 1;
iti.vol2(iti.vol2<1) = 1;
iti.vol3(iti.vol3<1) = 1;
iti.vol1(iti.vol1>2.5) = 3;
iti.vol2(iti.vol2>2.5) = 3;
iti.vol3(iti.vol3>2.5) = 3;
iti.vol1((iti.vol1~=1)&(iti.vol1~=3))=2;
iti.vol2((iti.vol2~=1)&(iti.vol2~=3))=2;
iti.vol3((iti.vol3~=1)&(iti.vol3~=3))=2;

stout.vol1(stout.vol1<1) = 1;
stout.vol2(stout.vol2<1) = 1;
stout.vol3(stout.vol3<1) = 1;
stout.vol1(stout.vol1>2.5) = 3;
stout.vol2(stout.vol2>2.5) = 3;
stout.vol3(stout.vol3>2.5) = 3;
stout.vol1((stout.vol1~=1)&(stout.vol1~=3))=2;
stout.vol2((stout.vol2~=1)&(stout.vol2~=3))=2;
stout.vol3((stout.vol3~=1)&(stout.vol3~=3))=2;

figure()
suptitle('Verification of the jitter values around a reversal (-3/+3 trials)')
subplot(3,3,1)
hist(iti.vol1,1:3)
ylabel('Occurence of the 3 iti jitters')
subplot(3,3,2)
hist(iti.vol2,1:3)
subplot(3,3,3)
hist(iti.vol3,1:3)
%
subplot(3,3,4)
hist(stout.vol1,1:3)
ylabel('Occurence of the 3 stout jitters')
subplot(3,3,5)
hist(stout.vol2,1:3)
subplot(3,3,6)
hist(stout.vol3,1:3)
%
subplot(3,3,7)
hist(give_combi(iti.vol1,stout.vol1),1:9)
ylabel('Occurence of the 9 iti-stout combinations')
xlabel('volatility level 1')
subplot(3,3,8)
hist(give_combi(iti.vol2,stout.vol2),1:9)
xlabel('volatility level 2')
subplot(3,3,9)
hist(give_combi(iti.vol3,stout.vol3),1:9)
xlabel('volatility level 3')

end