Check = cell(4,length(MiceOpt));

for midx = 1:length(MiceOpt)
    % compare
    tmpdir = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','UnitMatch.mat'));
    tmpMatchTbl = matfile(fullfile(tmpdir.folder,tmpdir.name));
    MatchTable = tmpMatchTbl.MatchTable;
    nclus = sqrt(height(MatchTable));
    MatchProbability = reshape(MatchTable.MatchProb,nclus,nclus);
    RankScore = reshape(MatchTable.RankScore,nclus,nclus);
    RankThreshold = reshape(MatchTable.SigFingerprintR,nclus,nclus);
    UMparam = tmpMatchTbl.UMparam;
    UniqueIDConversion = tmpMatchTbl.UniqueIDConversion;
    if UMparam.GoodUnitsOnly
        GoodId = logical(UniqueIDConversion.GoodID);
    else
        GoodId = true(1,length(UniqueIDConversion.GoodID));
    end
    UniqueID = UniqueIDConversion.UniqueID(GoodId);
    OriID = UniqueIDConversion.OriginalClusID(GoodId);
    OriIDAll = UniqueIDConversion.OriginalClusID;
    recses = UniqueIDConversion.recsesAll(GoodId);
    recsesall = UniqueIDConversion.recsesAll;

    % Individual scores:
    Scores2Inclde = UMparam.Scores2Include;
    % Add functional scores
    Scores2Inclde = {Scores2Inclde{:} 'ACGCorr','FRDiff','FingerprintCor'};
    for scid = 1:length(Scores2Inclde)
        eval([Scores2Inclde{scid} ' = reshape(MatchTable.' Scores2Inclde{scid} ',nclus,nclus);'])
    end

    for fid = 1:2
        % Now load manual curation
        if fid == 1
            tmpmanual = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','BlindFigures_WOFunct','BlindTable.mat'));
            tmpmanual = load(fullfile(tmpmanual.folder,tmpmanual.name));
            manualscore = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','BlindFigures_WOFunct','manualCuration*.mat'));
            if isempty(manualscore)
                continue
            end
            ScorerNames = cell(length(manualscore),2);
            AddName = 'WO';

            % Pairs of Interest
            tbl = tmpmanual.tbl;

            % UM results
            Pair1 = cell2mat(arrayfun(@(X) find(ismember(OriID,tbl.ClusID1(X)) & ismember(recses,tbl.RecID1(X))),1:height(tbl),'Uni',0));
            Pair2 = cell2mat(arrayfun(@(X) find(ismember(OriID,tbl.ClusID2(X)) & ismember(recses,tbl.RecID2(X))),1:height(tbl),'Uni',0));
            Pairs = arrayfun(@(X) [Pair1(X) Pair2(X)],1:length(Pair1),'Uni',0);
            MatchProb = cell2mat(cellfun(@(X) MatchProbability(X(1),X(2)),Pairs,'Uni',0))';
            Rank = cell2mat(cellfun(@(X) RankScore(X(1),X(2)),Pairs,'Uni',0))';
            RankThreshold = cell2mat(cellfun(@(X) RankThreshold(X(1),X(2)),Pairs,'Uni',0))';
            WithinSameSession = cell2mat(cellfun(@(X) recses(X(1)) == recses(X(2)),Pairs,'Uni',0));

            PyKS = cell2mat(cellfun(@(X) OriID(X(1)) == OriID(X(2)),Pairs,'Uni',0));

            % Individual scores
            for scid = 1:length(Scores2Inclde)
                eval([Scores2Inclde{scid} ' = cell2mat(cellfun(@(X) ' Scores2Inclde{scid} '(X(1),X(2)),Pairs,''Uni'',0));'])
            end

            Order =  1:height(tbl);
        else
            tmpmanual = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','BlindFigures','BlindTable.mat'));
            tmpmanual = load(fullfile(tmpmanual.folder,tmpmanual.name));
            tbl2 = tmpmanual.tbl;
            % Order?
            Order = cell2mat(arrayfun(@(X) find(tbl2.ClusID1 == tbl.ClusID1(X) & tbl2.ClusID2 == tbl.ClusID2(X) & tbl2.RecID1 == tbl.RecID1(X) & tbl2.RecID2 == tbl.RecID2(X),1,'first'),1:height(tbl),'Uni',0));
            ReOrderedTbl = tbl2(Order,:);
            manualscore = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','BlindFigures','manualCuration*.mat'));
            if isempty(manualscore)
                continue
            end
            AddName = 'With';

        end
        for id = 1:length(manualscore)

            tmpman = load(fullfile(manualscore(id).folder,manualscore(id).name));
            if id==1 && fid==1
                Manual = tmpman.match(Order);
            else
                Manual = cat(1,Manual,tmpman.match(Order));
            end
            tmpname = strsplit(manualscore(id).name,'manualCuration_');
            tmpname = strsplit(tmpname{2},'.mat');

            ScorerNames{id,fid} = [tmpname{1} AddName];
        end
    end
    ScorerNames = ScorerNames(~cellfun(@isempty,ScorerNames));

    %% How well do the different scoring methods correlate?
    AvgMan = nanmean(Manual,1);

    disp('User scoring correlations:')
    AllScoringMethods = cat(2,Manual',AvgMan',MatchProb,PyKS',Rank==1,RankThreshold);
    % Normalize between 0 and 1
    AllScoringMethods = (AllScoringMethods - nanmin(AllScoringMethods,[],1))./(nanmax(AllScoringMethods,[],1)-nanmin(AllScoringMethods,[],1));
    tmpcor = corr(AllScoringMethods);
    AllScorerNames = {ScorerNames{:},'AvgScorer','UM','KS','Rank1','RankTr'};
    if midx == 1
        corrfig = figure('name','Scorer correlations');
    else
        figure(corrfig)
    end
    subplot(ceil(sqrt(length(MiceOpt))),round(sqrt(length(MiceOpt))),midx)
    h=imagesc(tmpcor,[0,1]);
    set(gca,'XTick',1:size(tmpcor,1),'XTickLabel',{ScorerNames{:},'AvgScorer','UM','KS','Rank1','RankTr'},'YTick',1:size(tmpcor,1),'YTickLabel',{ScorerNames{:},'AvgScorer', 'UM','KS','Rank1','RankTr'})
    colormap(flipud(gray))
    title(MiceOpt{midx})
    colorbar
    axis square
    makepretty


    %% Any particular parameter correlations?
    if midx == 1
        % Save out correlation of each method with each of the parameters:
        CorrParam = nan(length(Scores2Inclde),length(AllScorerNames),length(MiceOpt));
        corrperparamfig = figure('name','Corerlations with parameters');
    else
        figure(corrperparamfig);
    end
    for scid = 1:length(Scores2Inclde)
        eval(['Idx = ~isnan(' Scores2Inclde{scid} ');'])
        eval(['CorrParam(scid,:,midx) = corr(' Scores2Inclde{scid} '(Idx)'',AllScoringMethods(Idx,:));'])
    end

    % Invert correlation for FR Diff (as less different = more likely
    % match)
    CorrParam(ismember(Scores2Inclde,'FRDiff'),:,midx) = -1*CorrParam(ismember(Scores2Inclde,'FRDiff'),:,midx);

    subplot(ceil(sqrt(length(MiceOpt))),round(sqrt(length(MiceOpt))),midx)
    h=imagesc(CorrParam(:,:,midx)',[0,1]);
    set(gca,'YTick',1:length(AllScorerNames),'YTickLabel',AllScorerNames,'XTick',1:length(Scores2Inclde),'XTickLabel',Scores2Inclde)
    colormap(flipud(gray))
    title(MiceOpt{midx})
    colorbar
    axis square
    makepretty

    

    % If PyKS said match
    disp('Compared to Stitched Kilosort:')
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(PyKS==1)>0)./sum(PyKS==1)*1000)/10) '% of PyKS Matches to be a match'])
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(PyKS==1)<0)./sum(PyKS==1)*1000)/10) '% of PyKS Matches to be a non-match'])
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(PyKS==0)>0)./sum(PyKS==0)*1000)/10) '% of PyKS Non-matches to be a match'])
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(PyKS==0)<0)./sum(PyKS==0)*1000)/10) '% of PyKS Non-matches to be a non-match'])

    % If UM said Match
    disp('Compared to UnitMatch:')
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(MatchProb>0.5)>0)./sum(MatchProb>0.5)*1000)/10) '% of UnitMatch Matches to be a match'])
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(MatchProb>0.5)<0)./sum(MatchProb>0.5)*1000)/10) '% of UnitMatch Matches to be a non-match'])
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(MatchProb<0.5)>0)./sum(MatchProb<0.5)*1000)/10) '% of UnitMatch Non-matches to be a match'])
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(MatchProb<0.5)<0)./sum(MatchProb<0.5)*1000)/10) '% of UnitMatch Non-matches to be a non-match'])


    if midx == 1
        clasfig = figure('name','Classifier Comparisons');
    else
        figure(clasfig)
    end
    subplot(2,2,1)
    if midx == 1
        line([0 length(MiceOpt)+1],[1 1],'color',[0 0 0],'LineStyle','--')
    end
    hold on
    clear h
    cols = lines(size(Manual,1));
    for id = 1:size(Manual,1)
        h(id) = scatter(midx+0.05*id,sum(Manual(id,PyKS==1)==1)./sum(PyKS==1),30,cols(id,:),'filled');
    end
    h(id+1) = scatter(midx,sum(MatchProb(PyKS==1)>0.5)./sum(PyKS==1),50,[0 0 0],'filled');
    h(id+2) = scatter(midx,sum(RankThreshold(PyKS==1)==1)./sum(PyKS==1),50,[0.75 0.75 0.75],'filled');

    title('Detection Performance matches relative to Kilosort Stitched')
    ylim([0 1])
    makepretty

    subplot(2,2,2)
    if midx == 1
        line([0 length(MiceOpt)+1],[1 1],'color',[0 0 0],'LineStyle','--')
    end
    hold on
    clear h
    cols = lines(size(Manual,1));
    for id = 1:size(Manual,1)
        h(id) = scatter(midx+0.05*id,sum(Manual(id,PyKS==0)==-1)./sum(PyKS==0),30,cols(id,:),'filled');
    end
    h(id+1) = scatter(midx,sum(MatchProb(PyKS==0)<0.5)./sum(PyKS==0),50,[0 0 0],'filled');
    h(id+2) = scatter(midx,sum(RankThreshold(PyKS==0)==0)./sum(PyKS==0),50,[0.75 0.75 0.75],'filled');

    title('Detection Performance non-matches relative to Kilosort Stitched')
    ylim([0 1])
    if midx == length(MiceOpt)
        legend([h(:)],{ScorerNames{:} 'UnitMatch','RankThrs'})
    end
    makepretty

    subplot(2,2,3)
    if midx == 1
        line([0 length(MiceOpt)+1],[1 1],'color',[0 0 0],'LineStyle','--')
    end
    hold on
    clear h
    cols = lines(size(Manual,1));
    for id = 1:size(Manual,1)
        h(id) = scatter(midx+0.05*id,sum(Manual(id,MatchProb>0.5)==1)./sum(MatchProb>0.5),30,cols(id,:),'filled');
    end
    h(id+1) = scatter(midx,sum(PyKS(MatchProb>0.5)==1)./sum(MatchProb>0.5),50,[0.5 0.5 0.5],'filled');
    h(id+2) = scatter(midx,sum(RankThreshold(MatchProb>0.5)==1)./sum(MatchProb>0.5),50,[0.75 0.75 0.75],'filled');

    title('Detection Performance matches relative to UnitMatch')
    ylim([0 1])
    makepretty

    subplot(2,2,4)
    if midx == 1
        line([0 length(MiceOpt)+1],[1 1],'color',[0 0 0],'LineStyle','--')
    end
    hold on
    clear h
    cols = lines(size(Manual,1));
    for id = 1:size(Manual,1)
        h(id) = scatter(midx+0.05*id,sum(Manual(id,MatchProb<0.5)==-1)./sum(MatchProb<0.5),30,cols(id,:),'filled');
    end
    h(id+1) = scatter(midx,sum(PyKS(MatchProb<0.5)==0)./sum(MatchProb<0.5),50,[0.5 0.5 0.5],'filled');
    h(id+2) = scatter(midx,sum(RankThreshold(MatchProb<0.5)==0)./sum(MatchProb<0.5),50,[0.75 0.75 0.75],'filled');

    title('Detection Performance non-matches relative to UnitMatch')
    ylim([0 1])
    if midx == length(MiceOpt)
        legend([h(:)],{ScorerNames{:} 'Kilosort','RankThrs'})
    end
    makepretty


    if midx == 1
        clasfigMan = figure('name','Classifier Comparisons to Average Manual Score');
    else
        figure(clasfigMan)
    end
    subplot(1,3,1)
    if midx == 1
        line([0 length(MiceOpt)+1],[1 1],'color',[0 0 0],'LineStyle','--')
    end
    hold on
    clear h
    cols = lines(size(Manual,1));
    for id = 1:size(Manual,1)
        h(id) = scatter(midx+0.05*id,sum(Manual(id,AvgMan>0.5)==1)./sum(AvgMan>0.5),30,cols(id,:),'filled');
    end
    h(id+1) = scatter(midx,sum(MatchProb(AvgMan>0.5)>0.5)./sum(AvgMan>0.5),50,[0 0 0],'filled');
    h(id+2) = scatter(midx-0.1,sum(PyKS(AvgMan>0.5)==1)./sum(AvgMan>0.5),50,[0.5 0.5 0.5],'filled');
    h(id+3) = scatter(midx-0.05,sum(RankThreshold(AvgMan>0.5)==1)./sum(AvgMan>0.5),50,[0.75 0.75 0.75],'filled');

    ylabel('Defined as match (proportion)')
    title(['Average manual score match'])
    ylim([0 1])

    xlabel('Dataset')

    makepretty

    subplot(1,3,2)
    if midx == 1
        line([0 length(MiceOpt)+1],[1 1],'color',[0 0 0],'LineStyle','--')
    end
    hold on
    clear h
    cols = lines(size(Manual,1));
    for id = 1:size(Manual,1)
        h(id) = scatter(midx+0.05*id,sum(Manual(id,AvgMan<-0.5)==1)./sum(AvgMan<-0.5),30,cols(id,:),'filled');
    end
    h(id+1) = scatter(midx,sum(MatchProb(AvgMan<-0.5)>0.5)./sum(AvgMan<-0.5),50,[0 0 0],'filled');
    h(id+2) = scatter(midx-0.1,sum(PyKS(AvgMan<-0.5)==1)./sum(AvgMan<-0.5),50,[0.5 0.5 0.5],'filled');
    h(id+3) = scatter(midx-0.05,sum(RankThreshold(AvgMan<-0.5)==1)./sum(AvgMan<-0.5),50,[0.75 0.75 0.75],'filled');

    ylabel('Defined as match (proportion)')

    title('Average manual score non-match')
    ylim([0 1])

    xlabel('Dataset')

    makepretty

    subplot(1,3,3)
    if midx == 1
        line([0 length(MiceOpt)+1],[1 1],'color',[0 0 0],'LineStyle','--')
    end
    hold on
    clear h
    cols = lines(size(Manual,1));
    for id = 1:size(Manual,1)
        h(id) = scatter(midx+0.05*id,sum(Manual(id,abs(AvgMan)<=0.5)==1)./sum(abs(AvgMan)<=0.5),30,cols(id,:),'filled');
    end
    h(id+1) = scatter(midx,sum(MatchProb(abs(AvgMan)<=0.5)>0.5)./sum(abs(AvgMan)<=0.5),50,[0 0 0],'filled');
    h(id+2) = scatter(midx-0.1,sum(PyKS(abs(AvgMan)<=0.5)==1)./sum(abs(AvgMan)<=0.5),50,[0.5 0.5 0.5],'filled');
    h(id+3) = scatter(midx-0.05,sum(RankThreshold(abs(AvgMan)<=0.5)==1)./sum(abs(AvgMan)<=0.5),50,[0.75 0.75 0.75],'filled');

    ylabel('Defined as match (proportion)')
    if midx == length(MiceOpt)
        legend([h(:)],{ScorerNames{:} 'UnitMatch','Kilosort','RankThrs'})
    end

    title('Average manual score uncertain')
    ylim([0 1])

    xlabel('Dataset')
    makepretty


    if midx == 1
        figFract = figure('name','Classifier Comparisons to Average Manual Score');
    else
        figure(figFract)
    end
    hold on
    clear h
    h(1) = scatter(midx,sum(AvgMan>0.5)./length(AvgMan),60,[0 0.5 0],'filled');
    h(2) = scatter(midx,sum(AvgMan<-0.5)./length(AvgMan),60,[0.5 0 0],'filled');
    h(3) = scatter(midx,sum(abs(AvgMan)<=0.5)./length(AvgMan),60,[0 0 0.5],'filled');
    title('Fraction')
    if midx == length(MiceOpt)
        legend(h,{'Everyone says match','Everyone says no match','Uncertain'})
    end

    %% How about if we split it up for within day and across day 'pairs'

    if midx == 1
        corrfigWithj = figure('name','Scorer correlations within session');
    else
        figure(corrfigWithj)
    end
    subplot(length(MiceOpt),2,(midx-1)*2+1)
    tmpcor = corr(AllScoringMethods(WithinSameSession,:));

    h=imagesc(tmpcor,[0,1]);
    if midx==length(MiceOpt)
        set(gca,'XTick',1:size(tmpcor,1),'XTickLabel',AllScorerNames,'YTick',1:size(tmpcor,1),'YTickLabel',AllScorerNames)
    else
        set(gca,'XTick',[],'XTickLabel',[],'YTick',1:size(tmpcor,1),'YTickLabel',AllScorerNames)

    end
    colormap(flipud(gray))
    title([MiceOpt{midx} ' within days'])
    colorbar
    axis square
    makepretty

    subplot(length(MiceOpt),2,midx*2)
    tmpcor = corr(AllScoringMethods(~WithinSameSession,:));

    h=imagesc(tmpcor,[0,1]);
    if midx == length(MiceOpt)
        set(gca,'XTick',1:size(tmpcor,1),'XTickLabel',AllScorerNames,'YTick',[],'YTickLabel',[])
    else
        set(gca,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[])

    end
    colormap(flipud(gray))
    title([MiceOpt{midx} ' across days'])
    colorbar
    axis square
    makepretty

    %% Sorting
     if midx == 1
        SortFig = figure('name','Sorting');
    else
        figure(SortFig)
     end
     subplot(2,length(MiceOpt),midx)
    
     tmpscores = AllScoringMethods(WithinSameSession,:);
     tmp = tmpscores(:,ismember(AllScorerNames,{'AvgScorer','UM','KS'}));
     tmp(:,1) = 20*(tmp(:,1)+0.1);
     [~,sortidx] = sort(nansum(tmp,2));
     tmpscores = cat(2,tmpscores(sortidx,strcmp(AllScorerNames,'AvgScorer')),tmpscores(sortidx,~strcmp(AllScorerNames,'AvgScorer')));
     h = imagesc(tmpscores);

     colormap redblue
     if midx == 1
         ylabel('Within recording')
     end
     title(MiceOpt{midx})
     set(gca,'XTick',1:length(AllScorerNames),'XTickLabel',{'AvgSc',AllScorerNames{~strcmp(AllScorerNames,'AvgScorer')}})
     makepretty


     subplot(2,length(MiceOpt),length(MiceOpt)+midx)
    
     tmpscores = AllScoringMethods(~WithinSameSession,:);
     tmp = tmpscores(:,ismember(AllScorerNames,{'AvgScorer','UM','KS'}));
     tmp(:,1) = 20*(tmp(:,1)+0.1);
     [~,sortidx] = sort(nansum(tmp,2));
     tmpscores = cat(2,tmpscores(sortidx,strcmp(AllScorerNames,'AvgScorer')),tmpscores(sortidx,~strcmp(AllScorerNames,'AvgScorer')));
     h = imagesc(tmpscores);
     colormap redblue
     if midx == 1
         ylabel('Across recordings')
     end
          set(gca,'XTick',1:length(AllScorerNames),'XTickLabel',{'AvgSc',AllScorerNames{~strcmp(AllScorerNames,'AvgScorer')}})

     makepretty

     %% Show which blind IDs to look at
     disp('BlindID (with functional scores) for which users think it''s a match, but UM not')
     idx = find(AvgMan' >0.5 & (MatchProb < 0.5));
     ReOrderedTbl.BlindID(idx)
     Check{1,midx} =  ReOrderedTbl.BlindID(idx);

     disp('BlindID (with functional scores) for which users think it''s not a match, but UM not')
     idx = find(AvgMan' <-0.5 & (MatchProb > 0.5));
     ReOrderedTbl.BlindID(idx) 
     Check{2,midx} =  ReOrderedTbl.BlindID(idx)

     disp('BlindID (with functional scores) for which users think it''s a match, but KS not')
     idx = find(AvgMan' >0.5 & (PyKS' == 0));
     ReOrderedTbl.BlindID(idx)
     Check{3,midx} =  ReOrderedTbl.BlindID(idx);

     disp('BlindID (with functional scores) for which users think it''s not a match, but KS not')
     idx = find(AvgMan' <-0.5 & (PyKS' == 1));
     ReOrderedTbl.BlindID(idx) 
     Check{4,midx} =  ReOrderedTbl.BlindID(idx)

end
midx = 5
Check{1,midx}(ismember(Check{1,midx},Check{3,midx}))
Check{2,midx}(ismember(Check{2,midx},Check{4,midx}))

idx = find(ismember(ReOrderedTbl.BlindID,37));
pairid1 = find(ismember(OriID,ReOrderedTbl.ClusID1(idx)) & ismember(recses,ReOrderedTbl.RecID1(idx)));
pairid2 = find(ismember(OriID,ReOrderedTbl.ClusID2(idx)) & ismember(recses,ReOrderedTbl.RecID2(idx)));
Pairs = {[pairid1 pairid2]};
pairid1 = find(ismember(OriID,ReOrderedTbl.ClusID2(idx)) & ismember(recses,ReOrderedTbl.RecID2(idx)));
pairid2 = find(ismember(OriID,ReOrderedTbl.ClusID1(idx)) & ismember(recses,ReOrderedTbl.RecID1(idx)));
Pairs{2} = [pairid1 pairid2];
MatchProb = cell2mat(cellfun(@(X) MatchProbability(X(1),X(2)),Pairs,'Uni',0))'
Rank = cell2mat(cellfun(@(X) RankScore(X(1),X(2)),Pairs,'Uni',0))'
RankThr = cell2mat(cellfun(@(X) RankThreshold(X(1),X(2)),Pairs,'Uni',0))'
WithinSameSession = cell2mat(cellfun(@(X) recses(X(1)) == recses(X(2)),Pairs,'Uni',0))

PyKS = cell2mat(cellfun(@(X) OriID(X(1)) == OriID(X(2)),Pairs,'Uni',0))