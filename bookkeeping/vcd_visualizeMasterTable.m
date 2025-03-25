function fH = vcd_visualizeMasterTable(master_table,store_imgs)

% set up figure basics
makeprettyfigures;

%%
strCols = [3,11:14,16]; % columns with names/characters (we don't want to plot)
colsToPlot = 1:size(master_table,2);
colsToPlot = setdiff(colsToPlot,strCols);

fH = figure(99); 
for jj = 1:length(unique(master_table.stim_class))
    
    stimClassData = master_table(master_table.stim_class==jj,:);
    taskClassData = stimClassData.task_class_name;

    nrows = length(colsToPlot);
    zoomwindow = 500;
    for ii = 1:nrows
        colName = master_table.Properties.VariableNames{colsToPlot(ii)};
        
        clf; set(gcf,'Position',[1,1,2500,1280]);
        
        dataToPlot = table2array(stimClassData(:,colsToPlot(ii)))';
        nsubplots = 1:ceil(length(dataToPlot)/500);
        
        if all(iscell(dataToPlot))
            dataToPlot = cell2mat(dataToPlot);
        end
        
        % convert to stim conditions
        if strcmp(colName,'contrast') || strcmp(colName,'rdk_coherence')
            nrconds0 = unique(dataToPlot);
            dataToPlot = sum((dataToPlot'==nrconds0).*[1:length(nrconds0)],2)';
            fixTickFlag = true;
        else
            fixTickFlag = false;
        end
        
        if  ~isequalwithequalnans(dataToPlot,NaN(size(dataToPlot)))
            cmap = cmapturbo(length(unique(dataToPlot)));
            
            % Plot all trials
            subplot(length(nsubplots)+1,1,1) 
            imagesc(dataToPlot);
            box off; 
            title(sprintf('%s - %s', stimClassData.stim_class_name{1}, colName),'Interpreter', 'none', 'FontSize', 25)
            
            % Set colorbar
            nrconds = unique(dataToPlot);
            nrconds = nrconds(~isnan(nrconds));
            colormap(cmap); cb = colorbar;
            if numel(nrconds)>1
                cmin = min(nrconds); cmax = max(nrconds);
            else
                cmin = 0; cmax = 1;
            end
            if fixTickFlag
                cb.Ticks = nrconds;
                cb.TickLabels = cellfun(@num2str, num2cell(nrconds0), 'UniformOutput', false);
            else
                if length(nrconds)<10
                    cb.Ticks = nrconds;
                elseif length(nrconds)<15
                    cb.Ticks = nrconds(1:2:end);
                else
                    cb.Ticks = nrconds(1:round(length(nrconds))/4:end);
                end
            end
            set(gca,'CLim', [cmin, cmax], 'YTick',[], 'XTick', [1:500:length(dataToPlot),length(dataToPlot)])
            ylabel('all trials');
            

            
            for kk = nsubplots
                trialidx = (kk-1)*zoomwindow  + [1:zoomwindow]';
                if trialidx(end)>length(dataToPlot)
                    trialidx(find(trialidx==(length(dataToPlot)+1)):end) = [];
                end
                if trialidx(1)<length(dataToPlot)
                    subplot(length(nsubplots)+1,1,kk+1); cla
                    
                    imagesc(dataToPlot(trialidx));
                    box off
                    
                    ylabel(sprintf('rows %d-%d', trialidx(1),trialidx(end)))
                    colormap(cmap); cb = colorbar;
                    
                    if fixTickFlag
                        cb.Ticks = nrconds;
                        cb.TickLabels = cellfun(@num2str, num2cell(nrconds0), 'UniformOutput', false);
                    else
                        if length(nrconds)<6
                            cb.Ticks = nrconds;
                        elseif length(nrconds)<10
                            cb.Ticks = nrconds(1:2:end);
                        else
                            cb.Ticks = nrconds(1:round(length(nrconds))/4:end);
                        end
                    end
                    set(gca,'CLim', [cmin, cmax],'YTick',[],'XTick',[1,[1:3].*round(length(trialidx)/4),length(trialidx)],...
                        'XTickLabel',trialidx([1,[1:3].*round(length(trialidx)/4),end]))
                    if kk == nsubplots
                        xlabel('Unique stimulus/trial nr');
                    end
                end
            end
            if store_imgs
                saveFigsFolder = fullfile(vcd_rootPath,'figs');
                filename = sprintf('vcd_stimclass%02d_%s.png', jj, colName);
                print(fH,'-dpng','-r300',fullfile(saveFigsFolder,filename));
            end
        end
    end
end