function [fit_table, RMSE] = PredictedWell(err_thresh,Ypred,Y,Ywt,Const,ConstLines,Gcl,GeneLabels)


Ypredmeans=splitapply(@(x) mean(x,1), Ypred', Gcl); % Predicted abundances averaged over replicates
Ymeans=Y'; % Experimental abundances averaged over replicates (from Box 1)

Error=100*(Ypredmeans-Ymeans)./Ywt'; % Calculate predicted error as a percentage of wt
M=size(Error,2); % # of transcripts and proteins

% Calculate the RMSE for each transcript and protein over all of the lines for each construct
RMSE=zeros(length(Const),M);
for c=1:length(Const)
    TargConst=Const{c};
    TargConstLine_ind=find(cellfun(@(x) contains(x,TargConst),ConstLines)); %indices for each line for the coi

    RMSE(c,:)=sqrt(sum(Error(TargConstLine_ind,:).^2,1))/length(TargConstLine_ind); 
end

fit_table=zeros(size(RMSE));
fit_table(RMSE<err_thresh)=1; % Designate a 'good' fit of the predicted abundance if its RMSE is < err_thresh

fit_table=array2table(fit_table,'VariableNames',GeneLabels,'RowNames',Const);