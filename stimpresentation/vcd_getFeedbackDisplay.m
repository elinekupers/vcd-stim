function [txt, text_rect] = vcd_getFeedbackDisplay(params, rect, behavioral_results)

% file path to movational phrases
perfphrase_folder = fullfile(params.instrfolder,'performance_phrases');

block_accuracy    = round(behavioral_results.summary.pct_correct);
mean_accuracy     = mean(behavioral_results.summary.pct_correct,'omitnan');
phrase_of_the_run = vcd_getMotivationText(mean_accuracy, perfphrase_folder);
txt = sprintf('%s\n',phrase_of_the_run);

for nn = 1:length(behavioral_results.summary.crossing_nr)
    
    taskName = strrep(params.exp.crossingnames{behavioral_results.summary.crossing_nr(nn)},'-','_');
    d = dir(fullfile(params.instrfolder,'instruction_txt',sprintf('%02d_runvcdcore_%s.txt', behavioral_results.summary.crossing_nr(nn), taskName)));
    
    script     = d(1).name;
    task_instr = vcd_getInstructionText(params, script, rect);
      
    task_name = task_instr(1: (regexp(task_instr, 'task')+3));
    
    if length(task_name)>20
        txt0 = sprintf('Block %d: %s \t %d%% correct', ...
            nn, task_name, block_accuracy(nn));
    elseif length(task_name)>17
        txt0 = sprintf('Block %d: %s \t\t %d%% correct', ...
            nn, task_name, block_accuracy(nn));
    elseif length(task_name)>10
        txt0 = sprintf('Block %d: %s \t\t\t %d%% correct', ...
            nn, task_name, block_accuracy(nn));
    else
        txt0 = sprintf('Block %d: %s \t\t\t\t %d%% correct', ...
            nn, task_name, block_accuracy(nn));
    end
    
    if length(txt0) < 60
        txt1 = strcat(txt0,repmat(' ',1,50-length(txt0)));
    else
        txt1 = txt0;
    end
    
    if nn~=length(behavioral_results.summary.crossing_nr)
       txt1 = strcat(txt1,'\n'); 
    end
    txt = [txt txt1];
        
end

% Create text rect
text_rect = rect + [params.offsetpix(1) params.offsetpix(2) params.offsetpix(1) params.offsetpix(2)];