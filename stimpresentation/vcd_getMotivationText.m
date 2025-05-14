function [instrtext, txt_rect, params] = vcd_getMotivationText(params, percentCorrect)

% function [instrtext, txt_rect, params] = vcd_getMotivationText(params, percentCorrect)
%
% Inputs:
% <params>         VCD parameters 
% <percentCorrect> integer or decimal, percentage of correct responses range 0-100
%
% Outputs:
% <instrtext>      char of the text to be shown
% <txt_file>       window rect 
% <params>         (optional?) updated params with .motivationalquotes_folder added which is
%                  the path where the text files live to motivate participants after a run
%
% Takes the percent correct and randomly grabs and loads a text file. There
% are 10 'good' (80% +); 10 'medium' (60-80%); 10 'bad' (60% or below); 
% and 1 'bug' if percentCorrect is empty or does not exist.
%
% Note use of dir -> will not work with spaces in the path

% file path to phrases
textfolder = fullfile(vcd_rootPath,'workspaces','instructions','performance_phrases');

% Percent correct into keyword
if exist('percentCorrect', 'var')
    if percentCorrect >= 80
        valence = 'good';
    elseif (percentCorrect > 60) && (percentCorrect < 80) % && error when bug: Operands to the logical and (&&) and or (||) operators must be convertible to logical scalar values.
        valence = 'medium';
    elseif percentCorrect <= 60 
        valence = 'bad';
    elseif isempty(percentCorrect)
        valence = 'bug';
    end
else   
    valence = 'bug';
end

% Find subset
prompts = dir([textfolder, sprintf('%s*.txt', valence)]);

% Randomly index into a .txt and set as txt_file 
text_ix = randi(length(prompts));
txt_file = strcat(prompts(text_ix).folder, '/', prompts(text_ix).name);

% Open text file
fd = fopen(txt_file, 'rt');

% Get characters
instrtext = '';
tl = fgets(fd);
lcount = 0;
while lcount < 48
    instrtext = [instrtext tl]; %#ok<*AGROW>
    tl = fgets(fd);
    lcount = lcount + 1;
end
% Close text file
fclose(fd);

% Create window rect
txt_rect = rect + [params.offsetpix(1) params.offsetpix(2) params.offsetpix(1) params.offsetpix(2)];

return