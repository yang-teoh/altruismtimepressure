filename1 = which(mfilename('fullpath'))
homepath= fileparts(filename1)

timelimit = [1.5 10];
subj = [1:16 18:33 35 36 38:60];
% participant 17, 34 & 37 poor calibration
where = NaN*zeros(10000, 16);
eyes2 =  NaN*zeros(10002, 160*length(subj));
for c = 1:2
    Eyes = zeros(timelimit(c)*1000, 80*length(subj))*NaN;
        for s = 1:length(subj)
        subject = num2str(subj(s))

        load(fullfile(homepath, 'PrimaryStudy','SubjectData',subject,['Data.' subject '.wfix.choice.mat']));

        selectedTrials = cellfun(@(x)logical(isequalfp(x, timelimit(c))), Data.TimeLimit);
        tempData.rt = cell2mat(Data.ChoiceRT(selectedTrials));

        Fix = Data.Fix(selectedTrials);


        maxrt = timelimit(c);


        for f = 1:length(Fix)
           currenteye = Fix{f};
           nfix = max(currenteye(:,2)); %number of fixations
    %        currenteye(nfix+1, 1) =  currenteye(nfix, 1);
    %        currenteye(nfix+1, 2) =  nfix + 1;
    %        currenteye(nfix+1, 3) =  currenteye(nfix, 4) + 1;
    %        currenteye(nfix+1, 4) =  round(maxrt*1000);
    %        currenteye(nfix+1, 5) =  currenteye(nfix, 5);
    %        currenteye(nfix+1, 6) =  currenteye(nfix, 6);
           for ey = 1:(nfix)%+1)
               p = currenteye(ey, 3);
                if p == 0
                    p = 1;
                end
            q = currenteye(ey, 4);
            if q > timelimit(c)*1000
                q = timelimit(c)*1000;
            end
            Eyes(p:q,(s-1)*80 + f) = currenteye(ey, 6);
           end
            Eyes((q+1):(timelimit(c)*1000),(s-1)*80 + f) = 0;
            eyes2(2,(c-1)*80*length(subj)+(s-1)*80 + f) = subj(s);

        end
        end
    
   
eyes2(1, ((c-1)*80*length(subj)+1):(c*80*length(subj))) = c;
eyes2(3:(2+timelimit(c)*1000), ((c-1)*80*length(subj)+1):(c*80*length(subj))) = Eyes(1:end,:);

for t = 1:(timelimit(c)*1000)
    where(t,(c-1)*8+1) = sum(Eyes(t,:)==1);
    where(t,(c-1)*8+2) = sum(Eyes(t,:)==2);
    where(t,(c-1)*8+3) = sum(isnan(Eyes(t,:)));
    where(t,(c-1)*8+4) = sum(Eyes(t,:)==0);
    
    where(t,(c-1)*8+5) = where(t,(c-1)*8+1)/sum(where(t,((c-1)*8+1):((c-1)*8+3)));
    where(t,(c-1)*8+6) = where(t,(c-1)*8+2)/sum(where(t,((c-1)*8+1):((c-1)*8+3)));
    where(t,(c-1)*8+7) = where(t,(c-1)*8+3)/sum(where(t,((c-1)*8+1):((c-1)*8+3)));
    where(t,(c-1)*8+8) = where(t,(c-1)*8+4)/sum(where(t,((c-1)*8+1):((c-1)*8+4)));
    
end

end


p1 = find(where(:,8) > .5)
fin1 = p1(1) -1

eyes22 = eyes2;
% 
% if(exist('permutationsummary.mat'))
% permut = load('permutationsummary.mat')
% permut = permut.permut
% currentperm = length(permut) + 1
% else
%     permut = []
%     currentperm = 1
% end

permut = []
currentperm = 1
parpool(16)
for nperm = currentperm:1000
nperm
eyes6 = eyes2(1:fin1+2,:);

for s = subj
    x = randperm(3);
    eyes7 = eyes6(3:end, eyes6(2,:)==s);
    nothing = isnan(eyes7);
    self = eyes7==1;
    other = eyes7==2;
    eyes7(nothing)= x(1);
    eyes7(self)= x(2);
    eyes7(other)= x(3);
    eyes7(eyes7==3)= NaN;
    eyes6(3:end, eyes6(2,:)==s) = eyes7;

end
eyes22 = eyes6;
    
shortonly = eyes2(1,:) == 1;
eyes3 = eyes22(:,shortonly);
longonly = eyes2(1,:) == 2;
eyes4 = eyes22(:,longonly);

tt1 =[];
parfor k = 1:fin1
    x = ~isnan(eyes3(k+2,:)) & ~(eyes3(k+2,:)==0);
    condi = eyes3(1,x)';
    subje = eyes3(2,x)';
    dats =  eyes3(k+2,x)';
    tb1 = table(condi, subje, dats, 'VariableNames',{'cond', 'subj', 'fix'});
    tb1.fix = tb1.fix-1; 
    glm = fitglme(tb1, 'fix ~ 1 + (1|subj)','Distribution', 'binomial', 'Link', 'logit');
    tt1(k,:) =  [glm.Coefficients.tStat(1) glm.Coefficients.pValue(1)];
end

tt1(:,3) = tt1(:,2);
tt1(:,4) = tt1(:,3) < .01;
%s2 = tt1(:,4) & tt1(:,1) < 0;
%tt1(:,5) = NaN;
%tt1(s1 ,5) = .9;
%tt1(s2 ,5) = .1;


short01 = find(tt1(:,4)==0);
short01where = [short01' fin1] - [1 short01'];
xx = short01where > 1;
short01stop = [short01' fin1];
short01start = [1 short01'];
shortsigduration = short01where(xx);
shortsigstart = short01start(xx);
shortsigend = short01stop(xx);
candTs = zeros(1);
for ind = 1:length(shortsigduration);
    st = shortsigstart(ind);
    ed = shortsigend(ind)-1;
    candTs(ind) = sum(tt1(st:ed, 1));
end

[valTs posTs] = max(abs(candTs));
permut(nperm, 1) = valTs;
permut(nperm, 2) = candTs(posTs);


tt2 = zeros(fin1,2);
parfor k = 1:fin1
    x = ~isnan(eyes4(k+2,:)) & ~(eyes4(k+2,:)==0);
    condi = eyes4(1,x)';
    subje = eyes4(2,x)';
    dats =  eyes4(k+2,x)';
    tb2 = table(condi, subje, dats, 'VariableNames',{'cond', 'subj', 'fix'});
    tb2.fix = tb2.fix-1; 
    glm = fitglme(tb2, 'fix ~ 1 + (1|subj)','Distribution', 'binomial', 'Link', 'logit');
    tt2(k,:) =  [glm.Coefficients.tStat(1) glm.Coefficients.pValue(1)];
end


tt2(:,3) = tt2(:,2);
tt2(:,4) = tt2(:,3) < .01;



long01 = find(tt2(:,4)==0);
long01where = [long01' fin1] - [1 long01'];
yy = long01where > 1;
long01stop = [long01' fin1];
long01start = [1 long01'];
longsigduration = long01where(yy);
longsigstart = long01start(yy);
longsigend = long01stop(yy);
candTl = zeros(1);
for ind = 1:length(longsigduration)
    st = longsigstart(ind);
    ed = longsigend(ind)-1;
    candTl(ind) = sum(tt2(st:ed, 1));
end



[valTl posTl] = max(abs(candTl));
permut(nperm, 3) = valTl;
permut(nperm, 4) = candTl(posTl);

eyes22 = eyes2;
for s = subj
    subcolumns = eyes22(2,:) == s;
    currentblocks = eyes22(1, subcolumns);
    x = randperm(length(currentblocks));
    eyes22(1,subcolumns) = currentblocks(x);
end

tt = zeros(fin1,2);

parfor j = 1:fin1
    x = ~isnan(eyes22(j+2,:)) & ~(eyes22(j+2,:)==0);
    tb = table(eyes22(1,x)', eyes22(2,x)', eyes22(j+2,x)');
    tb.Properties.VariableNames = {'cond', 'subj', 'fix'};
    tb.fix = tb.fix-1;
    glm = fitglme(tb, 'fix ~ cond  + (1|subj)','Distribution', 'binomial', 'Link', 'logit');
    tt(j,:) =  [glm.Coefficients.tStat(2) glm.Coefficients.pValue(2)];
end

tt(:,3) = tt(:,2);
tt(:,4) = tt(:,3) < .01;

lvs01 = find(tt(:,4)==0);
lvs01where = [lvs01' fin1]- [1 lvs01'];
zz = lvs01where > 1;
lvs01stop = [lvs01' fin1];
lvs01start = [1 lvs01'];
lvssigduration = lvs01where(zz);
lvssigstart = lvs01start(zz);
lvssigend = lvs01stop(zz);
candTlvs = zeros(1);
for ind = 1:length(lvssigduration)
    st = lvssigstart(ind);
    ed = lvssigend(ind)-1;
    candTlvs(ind) = sum(tt(st:ed, 1));
end


[valTlvs posTlvs] = max(abs(candTlvs));
permut(nperm, 5) = valTlvs;
permut(nperm, 6) = candTlvs(posTlvs);
end
delete(gcp('nocreate'))
save('clustercorrectpermutationsummary.mat','permut')

%load('permutationsummary.mat')
% 
% save('AvgAOIFixVisualisegrp.mat')
% load('AvgAOIFixVisualisegrp.mat')
% 
% figure()
% %subplot(3,1,1)
% plot( [1:(10000)],where(:,(1-1)*8+5), 'DisplayName', 'Self AOI','LineWidth', 3)
% hold on
% plot( [1:(10000)],where(:,(1-1)*8+6), 'Color', [0.9100    0.4100    0.1700],'DisplayName', 'Other AOI','LineWidth', 3)
% hold on
% plot( [1:10000],where(:,(1-1)*8+7),'LineStyle','--', 'Color', [.5,.5,.5],'DisplayName', 'Neither','LineWidth', 3)
% hold on
% plot( [1:fin1], tt1(:,5), 'r*', 'DisplayName', 'SelfvOther', 'MarkerSize',2)
% hold on
% plot([fin1,fin1], [0,1], 'k', 'LineWidth', 3)
% set(gca,'FontSize',12)
% xlim([0 5000])
% ylim([0 1])
% xlabel('Time(ms)')
% ylabel('Proportion Fixations')
% legend('Self AOI', 'Other AOI', 'Neither')
% 
% subplot(3,1,2)
% plot( [1:(10000)],where(:,(2-1)*8+5), 'DisplayName', 'Self AOI','LineWidth', 3)
% hold on
% plot( [1:(10000)],where(:,(2-1)*8+6),'Color', [0.9100    0.4100    0.1700],'DisplayName', 'Other AOI','LineWidth', 3)
% hold on
% plot( [1:10000],where(:,(2-1)*8+7),'LineStyle','--', 'Color', [.5,.5,.5],'DisplayName', 'Neither','LineWidth', 3)
% hold on
% plot( [1:fin2], tt2(:,5), 'r*', 'DisplayName', 'SelfvOther', 'MarkerSize',2)
% hold on
% plot([fin2,fin2], [0,1], 'k', 'LineWidth', 2)
% set(gca,'FontSize',12)
% xlim([0 5000])
% ylim([0 1])
% xlabel('Time(ms)')
% ylabel('Proportion Fixations')
% legend('Self AOI', 'Other AOI', 'Neither')
% 
% sc = zeros(10000, 2)
% 
% sc(:,1) = where(:,5) - where(:,6)
% sc(:,2) = where(:,8+5) - where(:,8+6)
% 
% subplot(3,1,3)
% plot( [1:(10000)],sc(:,1), 'Color', [0.8600    0.3700    0.1700] ,'DisplayName', 'Short Self Bias', 'LineWidth', 3)
% hold on
% plot( [1:(10000)],sc(:,2), 'Color', [0.1700    0.3700    0.8100],'DisplayName', 'Long Self Bias', 'LineWidth', 3)
% hold on
% plot([1:fin1], tt(:,5), 'r*','DisplayName', 'Pairwise T', 'MarkerSize', 2)
% hold on
% plot([fin1,fin1], [-1,1], 'k', 'LineWidth', 3)
% set(gca,'FontSize',12)
% 
% xlim([0 5000])
% ylim([-1 1])
% xlabel('Time(ms)')
% ylabel('Fixation Bias to Self')
% legend('Short Self Bias', 'Long Self Bias')
