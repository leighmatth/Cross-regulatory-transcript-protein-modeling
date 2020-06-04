function []=NodeDegrees(EdgeTable,GeneLabels)

M=length(GeneLabels);

InDegree=zeros(1,M);
InDeg_fromTran=zeros(1,M);
InDeg_fromProt=zeros(1,M);
OutDegree=zeros(1,M);
OutDeg_toTran=zeros(1,M);
OutDeg_toProt=zeros(1,M);
for i=1:M
    InDegree(1,i)=sum(strcmp(EdgeTable.Target,GeneLabels{i}));
    in_inds=find(strcmp(EdgeTable.Target,GeneLabels{i}));
    in_fromtype=cellfun(@(x) x(1), EdgeTable.Source(in_inds),'UniformOutput',0);
    InDeg_fromTran(1,i)=sum(strcmp(in_fromtype,'t'));
    InDeg_fromProt(1,i)=sum(strcmp(in_fromtype,'p'));
    
    OutDegree(1,i)=sum(strcmp(EdgeTable.Source,GeneLabels{i}));
    out_inds=find(strcmp(EdgeTable.Source,GeneLabels{i}));
    out_totype=cellfun(@(x) x(1), EdgeTable.Target(out_inds),'UniformOutput',0);
    OutDeg_toTran(1,i)=sum(strcmp(out_totype,'t'));
    OutDeg_toProt(1,i)=sum(strcmp(out_totype,'p'));
end

stackData_t=zeros(20,2,2);
stackData_t(:,1,1)=InDeg_fromTran(1:20);
stackData_t(:,1,2)=InDeg_fromProt(1:20);
stackData_t(:,2,1)=OutDeg_toTran(1:20);
stackData_t(:,2,2)=OutDeg_toProt(1:20);

ht=plotBarStackGroups(stackData_t, GeneLabels(1:20));
colors=lines(size(ht,2));
colors=[colors; colors+.5*(ones(2,3)-colors)];
colors=mat2cell(colors,ones(size(colors,1),1),3);
set(ht,{'FaceColor'},colors)
title('In and Out Degree of Transcripts','fontsize',22)
ylabel('# of Edges','fontsize',18)
set(gca,'XTickLabelRotation',45, 'YLim',[0 10])
set(gcf,'Color','w')
legend('In - from transcript','In - from protein', 'Out - to transcript', 'Out - to protein','Location','Best','fontsize',14)

stackData_p=zeros(20,2,2);
stackData_p(:,1,1)=InDeg_fromTran(21:40);
stackData_p(:,1,2)=InDeg_fromProt(21:40);
stackData_p(:,2,1)=OutDeg_toTran(21:40);
stackData_p(:,2,2)=OutDeg_toProt(21:40);

hp=plotBarStackGroups(stackData_p, GeneLabels(21:40));
set(hp,{'FaceColor'},colors)
title('In and Out Degree of Proteins','fontsize',22)
ylabel('# of Edges','fontsize',18)
set(gca,'XTickLabelRotation',45,'YLim',[0 10])
set(gcf,'Color','w')
legend('In - from transcript','In - from protein', 'Out - to transcript', 'Out - to protein','Location','Best','fontsize',14)
