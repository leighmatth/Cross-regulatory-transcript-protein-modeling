function Ypred=Model_Prediction(B,mue,Xtarg,Xmask,targ_prot_flag, Ywt)%Xact,Xtarg,Ywt)

if size(Xmask)~=size(Xtarg)
    error('Xmask and Xtarg must have same dimensions.')
end

if nargin < 6
    if targ_prot_flag==0
        Ywt = [];
    else
        error('Ywt must be passed in if targ_prot_flag is 1')
    end
else
    Ywt=Ywt{:,:}';
end

[M,Nact]=size(Xtarg);

if targ_prot_flag==1
    % Set proteins of targeted genes to percentage of their wild-type abundance
    Ytranswtavg=Ywt(1:M/2,1);
    Yprotwtavg=Ywt(1+M/2:M,1);
    Yprotnewtarg=Yprotwtavg.*(Xtarg(1:M/2,:)./Ytranswtavg);
    
    Xtarg=[Xtarg(1:M/2,:); Yprotnewtarg]; % adjust targeted values to include the proteins
    Xmask=[Xmask(1:M/2,:);Xmask(1:M/2,:)];
end
% Allocate space for variable
Ypred=zeros(M,Nact);

for j=1:Nact
    
    Kappa=diag(~Xmask(:,j)); %for experiment j, diagonal matrix with 0 for targeted transcripts/proteins, 1 for others
    
    Ypred(:,j)=(eye(M)-Kappa*B)\(Kappa*mue+Xtarg(:,j));
end