function EdgeTable = CreateEdgeTable(B,cutoff,select_table,Y,X_mask,Ypred,Ywt,Gc,Const,Gcl,ConstLines)
GeneLabels=Ywt.Properties.VariableNames;

jinds=zeros(length(Const),1);
jinds(1,1)=find(Gcl==1,1,'first'); %wt

for c=2:length(Const)
TargConst=Const{c}; % Construct
TargConst_ind=find(cellfun(@(x) strcmp(x,TargConst),Const)); %indices of the construct

% TargConstLine_ind=find(cellfun(@(x) contains(x,TargConst),ConstLines)); %indices for each line for the construct

% Find indices for targeted gene(s)
TargGene_inds=find(Gc==TargConst_ind); % All experiments
it=find(X_mask{TargGene_inds(1),:}); 

% Find the line/replicate that had the greated knockdown in the targeted
% gene transcripts
Ytargsub=Y(TargGene_inds,it); % Targeted abundances of the construct 
Ytargsub{:,:}=Ytargsub{:,:}./Ywt{:,it};
[~,min_rep]=min(sum(Ytargsub{:,:},2));

jinds(c,1)=TargGene_inds(min_rep); % experiment index of greatest targeted knockdown for this construct
end


EdgeTable=table(); % Initialize EdgeTable variable
M=size(Ypred,1); % Number of transcripts/proteins

for el=1:M

nonsigexps=find(select_table{:,el}==0); % Experiments where element el does not meet selection criteria (e.g., diff exp, 'good' fit)

col_inds=find(B(el,:)); % Edges influencing element el

parts=zeros(length(jinds),M);
for j=1:length(jinds)
    parts(j,:)=B(el,:).*Ypred(:,jinds(j))'; % Contribution of each influencing component on element el
end

parts_diff=parts(2:end,:)-repmat(parts(1,:),size(parts,1)-1,1); % Change in amount contributed in each experiment compared to contribution to wt abundance 

const_hm=parts_diff(:,col_inds); % pull out just the influencing elements
const_hm(nonsigexps,:)=[]; % Remove experiments where selection critia (e.g., diff exp, 'good' fit) is not met

totalchange_sigexps=sum(const_hm,2); % total change from wt for each experiment
perc_edge=const_hm./totalchange_sigexps; % amount each influencing element contributes as a function of total change

[~,c]=find(abs(perc_edge)>=cutoff); % Find influencing elements that contribute >= cutoff value (e.g., 50% of total change amount)

% Create EdgeTable with columns 'Source', 'Target', and 'Sign'
source=GeneLabels(col_inds(unique(c)))';
target=repmat(GeneLabels(el),length(source),1);
weight=sign(B(el,col_inds(unique(c))))';

EdgeTable=[EdgeTable; table(source,target,weight,'VariableNames',{'Source','Target','Sign'})];

end