function [DblEdges_in, DblEdges_out]=DblEdge_Motif(EdgeTable,GeneLabels)

M=length(GeneLabels);

DblEdges_in=table();
for i=1:M
    in_inds=find(strcmp(EdgeTable.Target,GeneLabels{i}));
    fromtype=cellfun(@(x) x(2:end),EdgeTable.Source(in_inds),'UniformOutput',0);
    [C, ia, ic]=unique(fromtype,'stable');
    if length(C)<length(fromtype)
        dup_ind=setdiff(1:numel(fromtype),ia);
        for j=1:length(dup_ind)
            % check that weights are opposite
            dup_inds_j=in_inds(find(ic==ic(dup_ind(j))));
            if sum(EdgeTable.Sign(dup_inds_j))==0
                DblEdges_in=[DblEdges_in; EdgeTable(dup_inds_j,:)];
            end
        end    
    end
end

DblEdges_out=table();
for i=1:M
    out_inds=find(strcmp(EdgeTable.Source,GeneLabels{i}));
    totype=cellfun(@(x) x(2:end),EdgeTable.Target(out_inds),'UniformOutput',0);
    [C, ia, ic]=unique(totype,'stable');
    if length(C)<length(totype)
        dup_ind=setdiff(1:numel(totype),ia);
        for j=1:length(dup_ind)
            % check that weights are opposite
            dup_inds_j=out_inds(find(ic==ic(dup_ind(j))));
            if sum(EdgeTable.Sign(dup_inds_j))==0
                DblEdges_out=[DblEdges_out; EdgeTable(dup_inds_j,:)];
            end
        end
            
    end
    
end