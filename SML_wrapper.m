%% SML wrapper function
function [B,mue,ilambda_cv,rho_factor,sigma2] = SML_wrapper(Y,Xmask,params)

K = ~logical(Xmask); % not Xmask to create K

% Use cross-validation to determine the SML regularization parameter
[ilambda_cv] = cross_validation_SML_B(Y,K,params);

%Scale each transcript and protein so that their abundances fall between 0
%and 1
[Ysc, Asc] = scaleY(Y);

% Use cross-validation to determine ridge regularization parameter and
% estimate sigma2 for SML
[rho_factor, sigma2] = cross_validation_ridge_B(Ysc,K,params);

% Calculate ridge regression of Ysc using rho_factor from the cross_validation
[BRhat, mueRhat] = constrained_ridge_B(Ysc,K,rho_factor);

% Initialize variables for SML
W = 1./abs(BRhat); % Use BRhat to compute the regularization weights for the SML algorithm
Q = Inf; % Ignore discarding rules in SML algorithm
lambda_factor_prev = 1;
BLhat = zeros(size(Y,1));

% Iterate SML algorithm through ilambdas until ilambda_cv
lambda_factors = params.lambda_factors; %SML regularization parameters
for ilambda = 1:ilambda_cv
    
    [BLhat, mueLhat] = sparse_maximum_likelihood_B(W,BLhat,Ysc,K,Q,lambda_factors(ilambda),lambda_factor_prev,sigma2,params);    
    
    Q = Inf; % Ignore discarding rules in SML algorithm
    lambda_factor_prev = lambda_factors(ilambda);
end %ilambda

% Use constrainedML algorithm on just the identified edges.
% This solves for the weights of the identifed edges without the l1-norm penalty
SL = abs(sign(BLhat)); % Detected edges from SML algorithm
[BChat, mueChat] = constrained_ML_B(BLhat,SL,Ysc,K,sigma2,params);

% Transform BChat and mueChat to use with un-scaled Y values
B = diag(1./Asc)*BChat*diag(Asc);
mue = diag(1./Asc)*mueChat;
sigma2 = sigma2*(1./Asc);