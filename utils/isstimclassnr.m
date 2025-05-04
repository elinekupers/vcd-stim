function var_out = isstimclassnr(var_in)
% Function to determine if input is a stimulus class numbers, regardless of
%  the class of the input,
%
% function var_out = isstimclassnr(var_in)
%
% <var_in> is a vector or cell
% if <var_in>, is an integer between 1 and 5, return <true>.  otherwise, return <false>.

input_type = class(var_in);

f_stimclassnr   = @(x) all(isnumeric(x) & (x>=1 & x <=5));


if strcmp(input_type,'cell')
    var_in = cell2mat(var_in);
elseif strcmp(input_type,'char')    
    var_in = str2num(var_in);
end

var_out = f_stimclassnr(var_in);
    
end