function var_out = choose_isempty(var_in)
% Function to determine if input is empty, regardless of the class of the
% input, given that isempty only works for arrays, not cells.
%
% function f = choose_isempty(var_in)
%
% <var_in> is a vector or cell
% if <var_in>, is empty, return <true>.  otherwise, return <false>.

input_type = class(var_in);

f_empty_array   = @(x) eval(sprintf('isempty([%s])',num2str(x)));
f_empty_cell    = @(x) isempty(cellfun("isempty",x));


if strcmp(input_type,'cell')
  
    if f_empty_cell(var_in)
        var_out = true;
    elseif ~f_empty_cell(var_in) && isequal(var_in{:},'')
        var_out = true;
    else
        var_out = false;
    end
    
elseif strcmp(input_type,'double')
    
    if f_empty_array(var_in)
        var_out = true;
    else
        var_out = false;
    end

elseif strcmp(input_type,'string')
    if strcmp(var_in,'')
        var_out = true;
    elseif isempty(var_in)
        var_out = true;
    else
        var_out = false;
    end
    
elseif strcmp(input_type,'char')
    if strcmp(var_in,"")
        var_out = true;
    elseif isempty(var_in)
        var_out = true;
    else
        var_out = false;
    end
    
end
    
end
