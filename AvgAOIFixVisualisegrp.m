filename1 = which(mfilename('fullpath'))
homepath= fileparts(filename1)
% timelimit = [1.5 10];
% subj =[1:16 18:33 35 36 38:60];
% %%% participant 17 34 and 37 poor calibration
% where = NaN*zeros(10000, 16);
% eyes2 =  NaN*zeros(10002, 160*length(subj));
% for c = 1:2
%     Eyes = zeros(timelimit(c)*1000, 80*length(subj))*NaN;
%         for s = 1:length(subj)
%         subject = num2str(subj(s))
% 
%         load(fullfile(homepath, 'PrimaryStudy','SubjectData',subject,['Data.' subject '.wfix.choice.mat']));
% 
%         selectedTrials = cellfun(@(x)logical(isequalfp(x, timelimit(c))), Data.TimeLimit);
%         tempData.rt = cell2mat(Data.ChoiceRT(selectedTrials));
% 
%         Fix = Data.Fix(selectedTrials);
% 
% 
%         maxrt = timelimit(c);
% 
% 
%         for f = 1:length(Fix)
%            currenteye = Fix{f};
%            nfix = max(currenteye(:,2)); %number of fixations
%            currenteye(nfix+1, 1) =  currenteye(nfix, 1);
%            currenteye(nfix+1, 2) =  nfix + 1;
%            currenteye(nfix+1, 3) =  currenteye(nfix, 4) + 1;
%            currenteye(nfix+1, 4) =  round(maxrt*1000);
%            currenteye(nfix+1, 5) =  currenteye(nfix, 5);
%            currenteye(nfix+1, 6) =  currenteye(nfix, 6);
%            for ey = 1:(nfix)%+1)
%                p = currenteye(ey, 3);
%                 if p == 0
%                     p = 1;
%                 end
%             q = currenteye(ey, 4);
%             if q > timelimit(c)*1000
%                 q = timelimit(c)*1000;
%             end
%             Eyes(p:q,(s-1)*80 + f) = currenteye(ey, 6);
%            end
%             Eyes((q+1):(timelimit(c)*1000),(s-1)*80 + f) = 0;
%             eyes2(2,(c-1)*80*length(subj)+(s-1)*80 + f) = s;
% 
%         end
%         end
%     
%    
% eyes2(1, ((c-1)*80*length(subj)+1):(c*80*length(subj))) = c;
% eyes2(3:(2+timelimit(c)*1000), ((c-1)*80*length(subj)+1):(c*80*length(subj))) = Eyes(1:end,:);
% 
% for t = 1:(timelimit(c)*1000)
%     where(t,(c-1)*8+1) = sum(Eyes(t,:)==1);
%     where(t,(c-1)*8+2) = sum(Eyes(t,:)==2);
%     where(t,(c-1)*8+3) = sum(isnan(Eyes(t,:)));
%     where(t,(c-1)*8+4) = sum(Eyes(t,:)==0);
%     
%     where(t,(c-1)*8+5) = where(t,(c-1)*8+1)/sum(where(t,((c-1)*8+1):((c-1)*8+3)));
%     where(t,(c-1)*8+6) = where(t,(c-1)*8+2)/sum(where(t,((c-1)*8+1):((c-1)*8+3)));
%     where(t,(c-1)*8+7) = where(t,(c-1)*8+3)/sum(where(t,((c-1)*8+1):((c-1)*8+3)));
%     where(t,(c-1)*8+8) = where(t,(c-1)*8+4)/sum(where(t,((c-1)*8+1):((c-1)*8+4)));
%     
% end
% 
% end
% 
% 
% p1 = find(where(:,8) > .5)
% fin1 = p1(1) -1
% p2 = find(where(:,16) > .5)
% fin2 = p2(1) - 1
% 
% f = eyes2(3:end,:);
% f(f==0) = NaN;
% f(f==2) = 0;
% eyes2(3:end,:) = f;
% 
% 
% shortonly = eyes2(1,:) == 1;
% eyes3 = eyes2(:,shortonly);
% longonly = eyes2(1,:) == 2;
% eyes4 = eyes2(:,longonly);
% 
% tt1 = zeros(fin1,2);
% for k = 1:fin1
%     x = ~isnan(eyes3(k+2,:));
%     tb1 = table(eyes3(1,x)', eyes3(2,x)', eyes3(k+2,x)');
%     tb1.Properties.VariableNames = {'cond', 'subj', 'fix'};
%     glm = fitglme(tb1, 'fix ~ 1 + (1|subj)','Distribution', 'binomial', 'Link', 'logit');
%     tt1(k,1) =  glm.Coefficients.tStat(1);
%     tt1(k,2) =  glm.Coefficients.pValue(1);
% end
% 
% 
% tt2 = zeros(fin2,2);
% for l = 412:fin1
%     x = ~isnan(eyes4(l+2,:));
%     tb1 = table(eyes4(1,x)', eyes4(2,x)', eyes4(l+2,x)');
%     tb1.Properties.VariableNames = {'cond', 'subj', 'fix'};
%     glm = fitglme(tb1, 'fix ~ 1 + (1|subj)','Distribution', 'binomial', 'Link', 'logit');
%     tt2(l,1) =  glm.Coefficients.tStat(1);
%     tt2(l,2) =  glm.Coefficients.pValue(1);
% end
% 
% tt = zeros(fin1,2);
% 
% for j = 1:fin1
%     x = ~isnan(eyes2(j+2,:));
%     tb = table(eyes2(1,x)', eyes2(2,x)', eyes2(j+2,x)');
%     tb.Properties.VariableNames = {'cond', 'subj', 'fix'};
%     glm = fitglme(tb, 'fix ~ cond  + (1|subj)','Distribution', 'binomial', 'Link', 'logit');
%     tt(j,1) =  glm.Coefficients.tStat(2);
%     tt(j,2) =  glm.Coefficients.pValue(2);
% end
% 
% 
% save('AvgAOIFixVisualisegrpexcl.mat')
load('AvgAOIFixVisualisegrpexcl.mat')

tt1(:,3) = tt1(:,2)
tt1(:,4) = tt1(:,3) < .01
s1 = tt1(:,4) & tt1(:,1) > 0
s2 = tt1(:,4) & tt1(:,1) < 0
tt1(:,5) = NaN
tt1(s1 ,5) = .9
tt1(s2 ,5) = .1

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




tt2(:,3) = tt2(:,2)
tt2(:,4) = tt2(:,3) < .01
s1 = tt2(:,4) & tt2(:,1) > 0
s2 = tt2(:,4) & tt2(:,1) < 0
tt2(:,5) = NaN
tt2(s1 ,5) = .9
tt2(s2 ,5) = .1

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




tt(:,3) = tt(:,2)
tt(:,4) = tt(:,3) < .01
s1 = tt(:,4) & tt(:,1) > 0
s2 = tt(:,4) & tt(:,1) < 0
tt(:,5) = NaN
tt(s1 ,5) = -.75
tt(s2 ,5) = .75


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



load('clustercorrectpermutationsummary.mat')

survshort = sum(abs(candTs) > permut(:,1)) > 950
survlong = sum(abs(candTl) > permut(:,3)) > 950
survlvs = sum(abs(candTlvs) > permut(:,5)) > 950

shortsigstart(survshort)
shortsigend(survshort)

longsigstart(survlong)
longsigend(survlong)

lvssigstart(survlvs)
lvssigend(survlvs)

tt1(:,5) = NaN
listsurvshort = find(survshort)
for j = listsurvshort
    st = shortsigstart(j)
    ed = shortsigend(j)
    if candTs(j) > 0 
    tt1(st:ed ,5) = .75
    else 
    tt1(st:ed ,5) = -.75
    end
end

tt2(:,5) = NaN
listsurvlong = find(survlong)
for j = listsurvlong
    st = longsigstart(j)
    ed = longsigend(j)
    if candTl(j) > 0 
    tt2(st:ed ,5) = .6
    else 
    tt2(st:ed ,5) = -.6
    end
end

tt(:,5) = NaN
listsurvlvs = find(survlvs)
for j = listsurvlvs
    st = lvssigstart(j)
    ed = lvssigend(j)
    if candTlvs(j) > 0 
    tt(st:ed ,5) = -.6
    else 
    tt(st:ed ,5) = .6
    end
end



figure()
set(gcf, 'DefaultAxesFontName', 'Arial');
set(gcf, 'DefaultTextFontName', 'Arial');
subplot(3,1,1)
plot( [1:fin1],where(1:fin1,(1-1)*8+5),'Color', [59/255,    113/255 ,   86/255], 'DisplayName', 'Self AOI','LineWidth', 1.25)
hold on
plot( [1:fin1],where(1:fin1,(1-1)*8+6), 'Color',[255/255,    150/255 ,   0],'DisplayName', 'Other AOI','LineWidth', 1.25)
hold on
plot( [1:fin1],where(1:fin1,(1-1)*8+7),'LineStyle','--', 'Color', [.5,.5,.5],'DisplayName', 'Neither','LineWidth', 1.25)
hold on
plot( [1:fin1], tt1(1:fin1,5), 'Color', 'k', 'DisplayName', 'SelfvOther')
hold on
plot([shortsigstart + (shortsigend-shortsigstart)/2 - 9, shortsigstart + (shortsigend-shortsigstart)/2 + 9], [.81, .81], 'k*') 
hold on
plot([fin1,fin1], [0,1], 'k', 'LineWidth', 2)
hold on
plot([shortsigend,shortsigend], [.725,.775], 'k')
hold on
plot([shortsigstart,shortsigstart], [.725,.775], 'k')
set(gca,'FontSize',8, 'FontName', 'Arial')
xlim([0 900])
ylim([0 1])
xlabel('Time(ms)')
ylabel({'Mean proportion fixations','under high time pressure'})
legend('Self AOI', 'Other AOI', 'Neither')
set(text(-100, 1.15, 'a'), 'FontSize', 8, 'FontName', 'Arial', 'FontWeight','Bold')

subplot(3,1,2)
plot( [1:fin1],where(1:fin1,(2-1)*8+5), 'Color', [59/255   , 113/255  ,  86/255],'DisplayName', 'Self AOI','LineWidth', 1.25)
hold on
plot( [1:fin1],where(1:fin1,(2-1)*8+6),'Color', [255/255  ,  150/255  ,  0],'DisplayName', 'Other AOI','LineWidth', 1.25)
hold on
plot( [1:fin1],where(1:fin1,(2-1)*8+7),'LineStyle','--', 'Color', [.5,.5,.5],'DisplayName', 'Neither','LineWidth', 1.25)
hold on
plot([fin1,fin1], [0,1], 'k', 'LineWidth', 2)
set(gca,'FontSize',8, 'FontName', 'Arial')
xlim([0 900])
ylim([0 1])
xlabel('Time(ms)')
ylabel({'Mean proportion fixations','under low time pressure'})
legend('Self AOI', 'Other AOI', 'Neither')
set(text(-100, 1.15, 'b'), 'FontSize', 8, 'FontName', 'Arial', 'FontWeight','Bold')

sc = zeros(10000, 2)

sc(:,1) = where(:,5) - where(:,6)
sc(:,2) = where(:,8+5) - where(:,8+6)

subplot(3,1,3)
plot( [1:fin1],sc(1:fin1,1), 'Color', [182/255  ,9/255,  84/255],'DisplayName', 'Short Self Bias', 'LineWidth', 1.25)
hold on
plot( [1:fin1],sc(1:fin1,2), 'Color', [64/255   , 144/255  ,  183/255],'DisplayName', 'Long Self Bias', 'LineWidth', 1.25)
hold on
plot([1:fin1], tt(1:fin1,5),'Color', 'k','DisplayName', 'Pairwise T')
hold on
plot([fin1,fin1], [-1,1], 'k', 'LineWidth', 2)
hold on
plot([lvssigend(1),lvssigend(1)], [.55,.65], 'k')
hold on
plot([lvssigstart(1),lvssigstart(1)], [.55,.65], 'k')
hold on
plot([lvssigstart(1) + (lvssigend(1) -lvssigstart(1))/2 - 9, lvssigstart(1) + (lvssigend(1) -lvssigstart(1))/2 + 9], [.73, .73], 'k*') 
hold on
plot([lvssigend(2),lvssigend(2)], [-.55,-.65], 'k')
hold on
plot([lvssigstart(2),lvssigstart(2)], [-.55,-.65], 'k')
hold on
plot([lvssigstart(2) + (lvssigend(2) -lvssigstart(2))/2 - 9, lvssigstart(2) + (lvssigend(2) -lvssigstart(2))/2 + 9], [-.73, -.73], 'k*')
hold on
plot([lvssigend(4),lvssigend(4)], [.55,.65], 'k')
hold on
plot([lvssigstart(4),lvssigstart(4)], [.55,.65], 'k')
hold on
plot([lvssigstart(4) + (lvssigend(4) -lvssigstart(4))/2 - 9, lvssigstart(4) + (lvssigend(4) -lvssigstart(4))/2 + 9], [.73, .73], 'k*')
set(gca,'FontSize',8, 'FontName', 'Arial')

xlim([0 900])
ylim([-1 1])
xlabel('Time(ms)')
ylabel({'Mean fixation bias', 'towards self outcomes'})
legend('High Time Pressure', 'Low Time Pressure', 'Location', 'southeast')
set(text(-100, 1.3, 'c'), 'FontSize', 8, 'FontName', 'Arial', 'FontWeight','Bold')
set(gcf, 'Units', 'centimeters')
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0,0 16.5, 17.5], 'PaperSize', [16.5, 17.5])
print('ATPconteyetracking', '-dpdf', '-r600')


