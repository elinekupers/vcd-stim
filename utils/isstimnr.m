function var_out = isstimnr(var_in)
% Function to determine if input is a valid stimulus numbers, regardless of
%  the class of the input
%
% function var_out = isstimnr(var_in)
%
% <var_in> is a vector or cell
% if <var_in>, is an integer between 1 and 1550, return <true>.  otherwise, return <false>.

input_type = class(var_in);

f_stimnr = @(x) all(x >=1 & x <= 1550);

if strcmp(input_type,'cell')
    if  ischar(var_in{1})
        var_in = str2double(str2mat(var_in))';
    else
        var_in = cell2mat(var_in);
    end
elseif strcmp(input_type,'char')    
    var_in = str2num(var_in);
end

var_out = f_stimnr(var_in) && isnumeric(var_in);
    
end