function [s] = table2html(t)
%     M = t(:)
%     sOut = makeHtmlTable(M, T, rowNames, colHeaders, colors, strPrecision)
%     
    % Convert table of unknown types to table of strings
    for ii=1:size(t,2) % Iterate through columns - each column can have a its own type
        columnName = t.Properties.VariableNames(1); % Get current column name
        tmp = t{:,1};
        switch(class(tmp))
            case 'categorical'
                tmp = cellstr(char(tmp));
            case 'double'
                tmp = cellstr(strtrim(num2str(tmp)));
            case 'logical'
                tmp = cellstr(bool2str(tmp));
            case 'cell'
            otherwise
                error('Unable to convert table column %d to string.',ii);
        end
        t(:,1) = []; % Remove pre-convert column
        t = [t table(tmp,'VariableNames',columnName)]; % Add converted column
    end
    
    % Build html table
    s = {};
    s{end+1} = sprintf('<html><table border="1"  cellpadding="4" cellspacing="0">\n');
    
     s{end+1} = sprintf('<tr>');
     
           if ~isempty(t.Properties.RowNames)
            
            s{end+1} = sprintf('<td>');
            s{end+1} = '';
            
            s{end+1} = sprintf('</td>');
            
        end
        for jj=1:size(t,2)
            s{end+1} = sprintf('<th>');
            s{end+1} = sprintf('%s', strtrim(char(t.Properties.VariableNames{jj})));
            s{end+1} = sprintf('</th>');
        end
        s{end+1} = sprintf('</tr>');
        
    for ii=1:size(t,1)
        s{end+1} = sprintf('<tr>');
        
        if ~isempty(t.Properties.RowNames)
            
            s{end+1} = sprintf('<td>');
            s{end+1} = sprintf('%s', strtrim(char(t.Properties.RowNames{ii})));
            
            s{end+1} = sprintf('</td>');
            
        end
        for jj=1:size(t,2)
            s{end+1} = sprintf('<td>');
            s{end+1} = sprintf('%s', strtrim(char(t{ii,jj})));
            s{end+1} = sprintf('</td>');
        end
        s{end+1} = sprintf('</tr>');
    end
    s{end+1} = sprintf('</table></html>\n');
end

% Convert boolean to string
function s = bool2str(x)
    s = [];
    for ii=1:numel(x)
        if x(ii)
            s(ii,:) = 'true ';
        else
            s(ii,:) = 'false';
        end
    end
end
