%% Box 1: Load and Prepare Data
load LigninTranscriptsProteins.mat

GeneNames = Y_all.Properties.VariableNames;
Experiments = Y_all.Properties.RowNames;

% Group the experiments by their experimental lines
Exp_Lines = cellfun(@(x) x(6:end-2),Experiments,'UniformOutput',0);
[G_lines, ID_lines] = findgroups(Exp_Lines);

% Take the average of the replicates for each line
Ytable = splitapply(@(x) mean(x,1),Y_all{:,:},G_lines);
Ytable = array2table(Ytable,'VariableNames',GeneNames,'RowNames',ID_lines);

% Each replicate for a line is the same so taking the average of the replicates
% will give a mask for the averagehas the same
Xmask_table = splitapply(@(x) mean(x,1),X_all_mask{:,:},G_lines);
Xmask_table = array2table(Xmask_table,'VariableNames',GeneNames,'RowNames',ID_lines);

% Transpose Y and Xmask so that each experiment represents a column and each row a
% transcript or protein
Y = Ytable{:,:}'; % {:,:} turns the table into an array
Xmask = Xmask_table{:,:}';

%% Box 2: Cross-validation
% Set seed for random number generator
rng_seed = 123456;
rng(rng_seed);

% Set the cross-validation parameters
params.Kcv = 5; % # of cross-validation folds
params.rho_factors = 10.^(-6:0.1:1); %Ridge regularization parameters
params.lambda_factors = 10.^(0:-.1:-3); %SML regularization parameters
params.maxiter = 1000; % max iterations for SML algorithm
params.cv_its = 1; % # of iterations of the kcv-fold cross-validations

%% Box 3
[B,mue,ilambda_cv,rho_factor,sigma2] = SML_wrapper(Y,Xmask,params);

%% Box 4: Check that the bounds of the regularization parameters were sufficient to find minimums
sml_regparam = [ilambda_cv > 1, ilambda_cv < length(params.lambda_factors)]; % Check if SML regularization parameter is either the first or last input parameter, indicating a minimum was not found
ridge_regparam = [rho_factor > params.rho_factors(1), rho_factor < params.rho_factors(end)]; % Check if ridge regression regularization parameter is either the first or last input parameter, indicating a minimum was not found

display(sml_regparam) % display result to command window; 1 if TRUE, 0 if FALSE
display(ridge_regparam) % display result to command window; 1 if TRUE, 0 if FALSE

%% Box 5: View/Compare result
Inferred_Relationships = abs(sign(B));
num_relationships = sum(sum(Inferred_Relationships)); % Number of relationships detected
display(num_relationships)

% Plot heatmap
colormap = [1 0 0; % red for negative relationships
    1 1 1; % white for no relationship detected
    0 1 0]; % green for positive relationships

figure;
heatmap(GeneNames,GeneNames,sign(B),'Colormap',colormap,'ColorbarVisible','off');
title('Inferred Relationships')

%% Box 6: Model Prediction Variables
Xmask = X_all_mask{:,:}'; % Mask indicating the targeted transcripts for all of the experiments
Xtarg = Xmask.*Y_all{:,:}'; % Desired knockdown levels of targeted transcripts
Ywt = Ytable{"WT_CK",:}'; % Average wildtype abundances
targ_prot_flag = 1; % Consider proteins targeted as well as transcripts

%% Box 7: Model Prediction
% Call the function to predict the un-targeted transcript and protein
% abundances
Ypred = Model_Prediction(B,mue,Xtarg,Xmask,targ_prot_flag,Ywt);
Ypred_table = array2table(Ypred','VariableNames',GeneNames,'RowNames',Experiments); % Save in table format with Gene names as columns and Experiments as rows

%% Box 8: Visualization
gene_of_interest='C3H3'; % Gene you want to plot
knockdown_of_interest='i29'; % Construct/experiment you want to plot

Knockdown_Visualization(gene_of_interest,knockdown_of_interest,Y_all,X_all_mask,Ypred)

%% Box 9: Determine transcripts and proteins that are predicted 'well' for each construct
err_thresh = 25; %error threshfold for being considered a 'good' fit

[Gc, Const] = findgroups(cellfun(@(x) x(6:8),Experiments,'UniformOutput',false)); %Group experiments into their different constructs (targeted knockdowns) (e.g., a01)
[Gcl, ConstLines] = findgroups(cellfun(@(x) x(6:end-2),Experiments,'UniformOutput',false)); %Group replicates of the same experimental line (e.g., a01.01 for construct a01 line 1)

[fit_table, RMSE] = PredictedWell(err_thresh,Ypred,Y,Ywt,Const,ConstLines,Gcl,GeneNames);


%% Box 10: Find the inferred relationships that contribute at least 50% of change from wt
load SigDEtable.mat
cutoff=0.5; %percentage of change from wt to identify highly influencing relationships (e.g., cutoff=0.5 -> 50% of total change amount)
select_table=fit_table(2:end,:); % creates table with columns and row names

select_table{:,:}=fit_table{2:end,:}.*SigDEtable{:,:}; % Look at experiments that were BOTH predicted 'well' and differentially expressed for each transcript/protein element

% Create EdgeTable with highly influencing relationships (aka edges)
EdgeTable = CreateEdgeTable(B,cutoff,select_table,Y_all,X_all_mask,Ypred,Ywt,Gc,Const,Gcl,GeneNames);

num_edges_subset = size(EdgeTable,1); % number of edges considered highly influencing
display(num_edges_subset)


%% Box 11: In/Out Degrees
NodeDegrees(EdgeTable,GeneNames) % Plots In-degrees (# of edges influencing a node (i.e., transcript/protein)) and Out-degrees (# of edges leaving a node)

%% Box 12: Network Motifs
% Find the double edge motifs
[DblEdges_in, DblEdges_out]=DblEdge_Motif(EdgeTable,GeneNames);

display(DblEdges_in)
display(DblEdges_out)