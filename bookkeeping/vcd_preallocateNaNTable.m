function nan_table = vcd_preallocateNaNTable(nr_rows, nr_cols, table_example, col_width)

% copy condition order table headers and scrub content
sz = [nr_rows nr_cols];
varTypes = varfun(@class,table_example,'OutputFormat','cell');

nan_table = table('Size',sz,'VariableTypes',varTypes, 'VariableNames',table_example.Properties.VariableNames);
for vt = 1:length(varTypes)
    
    if ~exist('col_width','var') || isempty(col_width) 
        % extract col_width: sometimes we deal with single table column, containing 2 sub columns
        cw = size(table_example.(table_example.Properties.VariableNames{vt}),2);
    elseif numel(col_width) > 1 && isequal(length(col_width),nr_cols)
        cw = col_width(vt);
    else
        cw = col_width;
    end
        
    if strcmp(varTypes(vt),'double')
        nan_table.(table_example.Properties.VariableNames{vt}) = NaN(sz(1),cw);
        
    elseif strcmp(varTypes(vt),'cell')
        if isequalwithequalnans(table2cell(table_example(1,vt)), {NaN(1,cw)})
            nan_table.(table_example.Properties.VariableNames{vt}) = repmat({NaN(1,cw)},sz(1),1);
        else
            nan_table.(table_example.Properties.VariableNames{vt}) = cell(sz(1),cw);
        end
    end
    
end

return