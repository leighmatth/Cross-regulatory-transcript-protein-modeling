% Copyright (c) 2019, North Carolina State University. All rights reserved.
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation under version 2 of the License.
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

function []=Knockdown_Visualization(g,TargConst,Yact_table,Xtable,Ypred)

colors=[180, 180, 180; % Gray 1
    204, 0, 0]./255; % Wolfpack Red


Experiments=Yact_table.Properties.RowNames;
GeneLabels=Yact_table.Properties.VariableNames;
Yact=Yact_table{:,:}';
[M,~]=size(Yact);

%%% Group experiments into their constructs
[Gc, Const]=findgroups(cellfun(@(x) x(6:8),Experiments,'UniformOutput',false));
[Gcl, ConstLines]=findgroups(cellfun(@(x) x(6:end-2),Experiments,'UniformOutput',false));

%%

i=find(strcmp(GeneLabels,['t' g])); %index of transcript for goi
ip=i+M/2; %index of protein for goi
TargConst_ind=find(cellfun(@(x) strcmp(x,TargConst),Const)); %indices of the coi

TargConstLine_ind=find(cellfun(@(x) contains(x,TargConst),ConstLines)); %indices for each line for the coi

WT_ind=find(cellfun(@(x) contains(x,'WT_'),Experiments));
Ywt=Yact_table{WT_ind,1:M}';
Ywt_avg=mean(Ywt,2);


TargGene_inds=find(Gc==TargConst_ind); % index for target gene
it=find(Xtable{TargGene_inds(1),:});
TargGeneNames=GeneLabels(it);

for n=1:length(TargConstLine_ind)
    jinds{n}=find(Gcl==TargConstLine_ind(n)); % indices coi
end


Ymeans=zeros(2,length(TargConstLine_ind)); % mean transcript of interest
Ymeans_p=zeros(2,length(TargConstLine_ind)); % mean protein of interest
Ymeans_t=zeros(length(it),length(TargConstLine_ind)); % mean targeted transcript(s)
for n=1:length(TargConstLine_ind)
    if numel(jinds{n})==1
        
        Ymeans_t(:,n)=Yact(it,jinds{n});
        Ymeans(:,n)=[Yact(i,jinds{n}) Ypred(i,jinds{n})];
        Ymeans_p(:,n)=[Yact(ip,jinds{n}) Ypred(ip,jinds{n})];
        
    else
        
        Ymeans_t(:,n)=mean(Yact(it,jinds{n})');
        Ymeans(:,n)=mean([Yact(i,jinds{n})' Ypred(i,jinds{n})']);
        Ymeans_p(:,n)=mean([Yact(ip,jinds{n})' Ypred(ip,jinds{n})']);
    end
end

% Calculate WTs
Ymeans_twt=mean([Yact(i,WT_ind)' Ypred(i,WT_ind)']);
Ymeans_pwt=mean([Yact(ip,WT_ind)' Ypred(ip,WT_ind)']);

% Combine Y_t and Y and Y_t and Y_tp
Ymean_trans=100*[Ymeans_twt; Ymeans']./Ywt_avg(i);
Ymean_prot=100*[Ymeans_pwt; Ymeans_p']./Ywt_avg(ip);

% sort lines by targeted knockdown level
[~, sort_it]=sort(sum(Ymeans_t./Ywt_avg(it),1),'descend');
sort_it=[1 1+sort_it];
Ymean_trans=Ymean_trans(sort_it,:);
Ymean_prot=Ymean_prot(sort_it,:);

wtjinds=[WT_ind jinds];
wtjinds=wtjinds(sort_it);

% plot transcript estimates
figure;
hold on

b=bar(Ymean_trans,'LineWidth',1);

b(1).FaceColor=[.9 .9 .9];
b(2).FaceColor=[245 194 190]./255;
b(1).EdgeColor=colors(1,:);
b(2).EdgeColor=colors(2,:);

pause(1)

for n=1:size(Ymean_trans,1)
    plot(b(1).XData(n)+b(1).XOffset,100*Yact(i,wtjinds{n})./Ywt_avg(i),'ko','MarkerFaceColor',colors(1,:),'MarkerSize',8)
    plot(b(2).XData(n)+b(2).XOffset,100*Ypred(i,wtjinds{n})./Ywt_avg(i),'ko','MarkerFaceColor',colors(2,:),'MarkerSize',8)
end

wt_line=refline(0,100);
wt_line.LineStyle='--';
wt_line.LineWidth=0.75;
wt_line.Color='k';
ylabel('Abundance (% of WT)','FontSize',18)
set(gca,'XTick',1:length(TargConstLine_ind)+1,'XTickLabel',['WT'; ConstLines(TargConstLine_ind(sort_it(2:end)-1))],'FontSize',18)
set(gcf,'Color','w')
ylimits=ylim;
if ylimits(2)<120
    set(gca,'YLim',[0 120])
end
title(GeneLabels(i),'FontSize',22)
leg_text={'Measured','Predicted'};
legend(leg_text,'Location','Best')


% Plot Protein Estimates
figure;
hold on

b=bar(Ymean_prot,'LineWidth',1);

b(1).FaceColor=[.9 .9 .9];
b(2).FaceColor=[245 194 190]./255;
b(1).EdgeColor=colors(1,:);
b(2).EdgeColor=colors(2,:);

pause(1)

for n=1:size(Ymean_trans,1)
    plot(b(1).XData(n)+b(1).XOffset,100*Yact(ip,wtjinds{n})./Ywt_avg(ip),'ko','MarkerFaceColor',colors(1,:),'MarkerSize',8)
    plot(b(2).XData(n)+b(2).XOffset,100*Ypred(ip,wtjinds{n})./Ywt_avg(ip),'ko','MarkerFaceColor',colors(2,:),'MarkerSize',8)
end

wt_line=refline(0,100);
wt_line.LineStyle='--';
wt_line.LineWidth=0.75;
wt_line.Color='k';
ylabel('Abundance (% of WT)','FontSize',18)
set(gca,'XTick',1:length(TargConstLine_ind)+1,'XTickLabel',['WT'; ConstLines(TargConstLine_ind(sort_it(2:end)-1))],'FontSize',18)
set(gcf,'Color','w')
ylimits=ylim;
if ylimits(2)<120
    set(gca,'YLim',[0 120])
end
title(GeneLabels(ip),'FontSize',22)
leg_text={'Measured','Predicted'};
legend(leg_text,'Location','Best')


%% Plot Targeted Bar plots
colors_gray=gray(3*length(it)+2);

Ymeans_targ=100*[Ywt_avg(it)'; Ymeans_t']./Ywt_avg(it)';
Ymeans_targ=Ymeans_targ(sort_it,:);

% Plot Transcript Estimates
figure;
hold on
b=bar(Ymeans_targ,'LineWidth',1);
for nt=1:length(it)
    b(nt).FaceColor=colors_gray(2*length(it)+nt+1,:);
    b(nt).EdgeColor=colors_gray(3*nt+1,:);
end

pause(1)
for nt=1:length(it)
    for n=1:size(Ymean_trans,1)
        plot(b(nt).XData(n)+b(nt).XOffset,100*Yact(it(nt),wtjinds{n})./Ywt_avg(it(nt)),'ko','MarkerFaceColor',colors_gray(3*nt+1,:),'MarkerSize',8)
    end
end
wt_line=refline(0,100);
wt_line.LineStyle='--';
wt_line.LineWidth=0.75;
wt_line.Color='k';
ylabel('Abundance (% of WT)','FontSize',18)
set(gca,'XTick',1:length(TargConstLine_ind)+1,'XTickLabel',['WT';ConstLines(TargConstLine_ind(sort_it(2:end)-1))],'FontSize',18)
set(gcf,'Color','w')
ylimits=ylim;
if ylimits(2)<120
    set(gca,'YLim',[0 120])
end
title('Transcripts of Targeted Genes','FontSize',22)
legend(TargGeneNames,'Location','Best')

