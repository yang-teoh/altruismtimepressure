%homepath = manual input

filename1 = which(mfilename('fullpath'))
homepath= fileparts(filename1)

for s = [1:60]

sub = num2str(s)
load(fullfile(homepath, 'PrimaryStudy','SubjectData', sub,[ 'Data.' sub '.choice.mat']));
p = dir([homepath filesep 'PrimaryStudy' filesep 'FIX' filesep 'Output' filesep ['Untitled_' sub '_*.csv']]);
Eye = readtable(fullfile(homepath, 'PrimaryStudy', 'FIX','Output', p.name));

Eye.IP_RELATIME_START = Eye.IP_START_TIME - Eye.TRIAL_START_TIME;
Eye.IP_RELATIME_END = Eye.IP_END_TIME - Eye.TRIAL_START_TIME;
Eye.CURRENT_FIX_TRIALTIME_START = Eye.CURRENT_FIX_START - Eye.IP_RELATIME_START;
Eye.CURRENT_FIX_TRIALTIME_END = Eye.CURRENT_FIX_END - Eye.IP_RELATIME_START;

j  = strcmp(Eye.CURRENT_FIX_INTEREST_AREA_ID, '.');
q = length(Eye.CURRENT_FIX_INTEREST_AREA_ID);

repl = [1:q];
repl = repl(j);
for i = repl;
Eye.CURRENT_FIX_INTEREST_AREA_ID{i} = NaN;
end

j  = strcmp(Eye.CURRENT_FIX_INTEREST_AREA_ID, '1');
q = length(Eye.CURRENT_FIX_INTEREST_AREA_ID);

repl = [1:q];
repl = repl(j);
for i = repl;
Eye.CURRENT_FIX_INTEREST_AREA_ID{i} = 1;;
end

j  = strcmp(Eye.CURRENT_FIX_INTEREST_AREA_ID, '2');
q = length(Eye.CURRENT_FIX_INTEREST_AREA_ID);

repl = [1:q];
repl = repl(j);
for i = repl;
Eye.CURRENT_FIX_INTEREST_AREA_ID{i} = 2;
end

j  = strcmp(Eye.CURRENT_FIX_BLINK_AROUND, 'NONE');
j = ~j;
q = length(Eye.CURRENT_FIX_BLINK_AROUND);
repl = [1:q];
for i = repl;
    if (any(i == repl(j)))
    Eye.CURRENT_FIX_BLINK_AROUND{i} = 1;
    else
    Eye.CURRENT_FIX_BLINK_AROUND{i} = 0;
    end
end

Eye.CURRENT_FIX_INTEREST_AREA_ID = cell2mat(Eye.CURRENT_FIX_INTEREST_AREA_ID);
Eye.CURRENT_FIX_BLINK_AROUND = cell2mat(Eye.CURRENT_FIX_BLINK_AROUND);

Fixate = zeros(length(Eye.CURRENT_FIX_INTEREST_AREA_ID), 7);
Fixate(:,1) = Eye.TRIAL_INDEX;
Fixate(:,2) = Eye.CURRENT_FIX_INDEX;
Fixate(:,3) = Eye.CURRENT_FIX_TRIALTIME_START;
Fixate(:,4) = Eye.CURRENT_FIX_TRIALTIME_END;
Fixate(:,5) = Eye.IP_RELATIME_END - Eye.IP_RELATIME_START;
Fixate(:,6) = Eye.CURRENT_FIX_INTEREST_AREA_ID;
Fixate(:,7) = Eye.CURRENT_FIX_BLINK_AROUND;
% Fixate = Fixate(~Fixate(:,7),:)
left = Fixate(:,6)==1;
right = Fixate(:,6)==2;
if (mod(s-1,4) + 1 < 3) ;
    Fixate(left,6) = 1;
    Fixate(right,6) = 2;
else
    Fixate(left,6) = 2;
    Fixate(right,6) = 1;
end

Data.Fix = cell(1, 160);

for i = 1:160
    ind = Fixate(:,1) == i;
    Data.Fix{i} = Fixate(ind,:);
end

save(fullfile(homepath, 'PrimaryStudy', 'SubjectData', sub,[ 'Data.' sub '.wfix.choice.mat']), 'Data')

end
