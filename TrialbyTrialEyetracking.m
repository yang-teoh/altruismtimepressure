dataframe = zeros(60*160, 10)
filename1 = which(mfilename('fullpath'))
homepath= fileparts(filename1)

for s = [1:60]
    subject = num2str(s)
    load(fullfile(homepath, 'PrimaryStudy','SubjectData',subject,['Data.' subject '.wfix.choice.mat']))

    NonResp = cellfun(@(x)strcmp(x,'NULL'),Data.Resp);
    Data.Resp(NonResp) = {NaN}; 
    Data.Resp = +(cell2mat(Data.Resp) == 2);
    Data.Resp(NonResp) = NaN;

    Data.SelfProposal = cell2mat(Data.SelfProposal);
    Data.OtherProposal = cell2mat(Data.OtherProposal);
    Data.TimeLimit = cell2mat(Data.TimeLimit);
    Data.ChoiceRT = cell2mat(Data.ChoiceRT);

    for i = 1:length(Data.Resp)
     dataframe((s-1)*160+i, 1) = s; %Subject ID
     dataframe((s-1)*160+i, 2) = i; %trial number
     dataframe((s-1)*160+i, 3) = Data.TimeLimit(i); %condition (time limit)
     dataframe((s-1)*160+i, 4) = Data.SelfProposal(i); %self
     dataframe((s-1)*160+i, 5) = Data.OtherProposal(i); %other
     dataframe((s-1)*160+i, 6) = Data.SelfProposal(i) - Data.OtherProposal(i); %ineq
     dataframe((s-1)*160+i, 7) = Data.Resp(i); %accept/reject
     dataframe((s-1)*160+i, 8) = Data.ChoiceRT(i); %rt
     
     temp4 =  cell2mat(Data.Fix(i)); %get fixation 
     if isempty(temp4) 
        dataframe((s-1)*160+i, 11) = 0; %no fixations early or late
        dataframe((s-1)*160+i, 12) = 0;
        dataframe((s-1)*160+i, 13) = 0;
        dataframe((s-1)*160+i, 14) = 0;
     else
         while (temp4(1,4) < 1) %catching bugs.
             temp4 = temp4(2:end,:);
         end
         if isempty(temp4)
             dataframe((s-1)*160+i, 11) = 0; 
             dataframe((s-1)*160+i, 12) = 0; 
             dataframe((s-1)*160+i, 13) = 0;
             dataframe((s-1)*160+i, 14) = 0; 
         else
             temp4(1,3) = 1
             if isnan(Data.ChoiceRT(i))
                veclen = Data.TimeLimit(i)*1000; %length of vector of gaze positions
             else
                veclen = ceil(Data.ChoiceRT(i)*1000); 
             end
             vec = zeros(1,veclen)*NaN; 
             %Filling in the vector with the AOIs 
            for j = 1:length(temp4(:,1))
                 st = temp4(j, 3);
                 ed = temp4(j, 4);
                 if ed > Data.TimeLimit(i)*1000;
                     ed = Data.TimeLimit(i)*1000;
                 end
                 vec(st:ed) = temp4(j,6);
            end
            if veclen < 286 %termination of the trial before this early period
                dataframe((s-1)*160+i, 11) = sum(vec(1:end) == 1); %early gaze proportion is calculated up to the point of termination
                dataframe((s-1)*160+i, 12) = sum(vec(1:end) == 2);
                dataframe((s-1)*160+i, 13) = 0;
                dataframe((s-1)*160+i, 14) = 0;
            else
                
                dataframe((s-1)*160+i, 11) = sum(vec(1:286) == 1);
                dataframe((s-1)*160+i, 12) = sum(vec(1:286) == 2);
                dataframe((s-1)*160+i, 13) = sum(vec(287:end) == 1);
                dataframe((s-1)*160+i, 14) = sum(vec(287:end) == 2);
            end
         end
        end
     
     if (Data.Resp(i)==1) 
         dataframe((s-1)*160+i, 15) = Data.OtherProposal(i)-50; %outcome for the other person
     elseif (Data.Resp(i)==0) 
         dataframe((s-1)*160+i, 15) = 50-Data.OtherProposal(i);
     end
     
     end
end
x = dataframe(:, 1) == 0
dataframe(x, :) = []
T = array2table(dataframe,'VariableNames', {'subj','trial','tp','self','other','fair','resp','rt','ffix','fdur','gazeself','gazeother', 'gazelateself','gazelateother','given'})
% csvwrite('TPEbehaveyetracking286.csv', dataframe)
writetable(T, 'PrimaryStudyEyetracking.csv')