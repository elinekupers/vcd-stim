function [instrtext, txt_rect] = vcd_getInstructionText(params, txt_file, rect)

% Open text file
fd = fopen(txt_file, 'rt');

% get characters
instrtext = '';
tl = fgets(fd);
lcount = 0;
while lcount < 48
    instrtext = [instrtext tl]; %#ok<*AGROW>
    tl = fgets(fd);
    lcount = lcount + 1;
end
% close text file
fclose(fd);

% create window rect
txt_rect = rect + [params.offsetpix(1) params.offsetpix(2) params.offsetpix(1) params.offsetpix(2)];

return