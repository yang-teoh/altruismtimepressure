dataframe = zeros(60*160, 10)
filename1 = which(mfilename('fullpath'))
homepath= fileparts(filename1)

for s = [1:60]
    subject = num2str(s)
   load(fullfile(homepath, 'PrimaryStudy','SubjectData',subject,['Data.' subject '.choice.mat']))

    NonResp = cellfun(@(x)strcmp(x,'NULL'),Data.Resp);
    Data.Resp(NonResp) = {NaN}; 
    Data.Resp = +(cell2mat(Data.Resp) == 2);
    Data.Resp(NonResp) = NaN;

    Data.SelfProposal = cell2mat(Data.SelfProposal);
    Data.OtherProposal = cell2mat(Data.OtherProposal);
    Data.TimeLimit = cell2mat(Data.TimeLimit);
    Data.ChoiceRT = cell2mat(Data.ChoiceRT);

    for i = 1:length(Data.Resp)
     dataframe((s-1)*160+i, 1) = s;
     dataframe((s-1)*160+i, 2) = 1;
     dataframe((s-1)*160+i, 3) = Data.TimeLimit(i);
     dataframe((s-1)*160+i, 4) = Data.SelfProposal(i);
     dataframe((s-1)*160+i, 5) = Data.OtherProposal(i);
     dataframe((s-1)*160+i, 6) = Data.SelfProposal(i) - Data.OtherProposal(i);
     dataframe((s-1)*160+i, 7) = Data.Resp(i);
     dataframe((s-1)*160+i, 8) = Data.ChoiceRT(i);
     if (Data.Resp(i)==1) 
         dataframe((s-1)*160+i, 9) = Data.OtherProposal(i)-50;
     elseif (Data.Resp(i)==0) 
         dataframe((s-1)*160+i, 9) = 50-Data.OtherProposal(i);
     end
    end
end

x = dataframe(:, 1) == 0
dataframe(x, :) = []
T = array2table(dataframe,'VariableNames',{'subj','exp','tp','self','other','fair','resp','rt','given','gen'})

writetable(T, 'PrimaryStudyBehav.csv')