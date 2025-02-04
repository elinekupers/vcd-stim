total_task = cell(1,length(trials.ns.tasks));


for ii = 1:length(trials.ns.tasks)
    total_n = [];
    for jj = 1:length(trials.ns.tasks(ii).miniblock)
        [a,b] = histcounts([trials.ns.tasks(ii).miniblock(jj).trial(:).super_cat],'BinEdges',[1:6]);
        total_n = cat(1,total_n,a);
    end
    total_task{ii} = total_n;
end


unique(sum(total_n,2))

%%
figure(1); clf;
for n = 1:length(trials.ns.tasks)
    subplot(1,length(trials.ns.tasks),n) 
    imagesc(total_task{n}); colorbar;
    set(gca,'CLim',[0 5])
    title(trials.ns.tasks(n).name)
end


