function nan_table = vcd_preallocateNaNTable(nr_rows, nr_cols, table_example, col_width)

% copy condition order table headers and scrub content
sz = [nr_rows nr_cols];
varTypes = varfun(@class,table_example,'OutputFormat','cell');
varNames = table_example.Properties.VariableNames;
% check matlab version for backwards compatibility
% if this version isn't from the 2000s, we need to do make the table the annoying way
if verLessThan('matlab', '9.6')
    
    nan_table = table();
    for ii = 1:length(varTypes)
        if strcmp(varTypes{ii},'cell')
            tmp = table(mat2cell(NaN(sz(1),1), ones(sz(1),1)), 'VariableNames', varNames(ii));
        elseif strcmp(varTypes{ii},'logical')
            tmp = table(false(sz(1),1), 'VariableNames', varNames(ii));
        elseif strcmp(varTypes{ii},'double')
            tmp = table(NaN(sz(1),1), 'VariableNames', varNames(ii));
        else
            error('[%s]: Not sure what to do with class type of %s.',mfilename,varTypes{ii})
        end
        nan_table = cat(2, nan_table, tmp);
    end
    
else
    % Otherwise we start the way we intend to, creating a "default" table 
    % (zeros for all rows) and update class type below. 
    nan_table = table('Size',sz,'VariableTypes',varTypes, 'VariableNames',varNames);
end

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
        
    elseif strcmp(varTypes(vt),'logical')
        nan_table.(table_example.Properties.VariableNames{vt}) = false(sz(1),cw);
        
    elseif strcmp(varTypes(vt),'cell')
        if isequalwithequalnans(table2cell(table_example(1,vt)), {NaN(1,cw)})
            nan_table.(table_example.Properties.VariableNames{vt}) = repmat({NaN(1,cw)},sz(1),1);
        else
            nan_table.(table_example.Properties.VariableNames{vt}) = cell(sz(1),cw);
        end
    else
        error('[%s]: Not sure what to do with class type of %s.',mfilename,varTypes{vt})
    end
    
end

return