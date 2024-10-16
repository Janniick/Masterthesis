function results_row = addNestedFieldsToResultsRow(structure, parentName, results_row)
    % Get all field names in the current structure
    fields = fieldnames(structure);
    
    for i = 1:length(fields)
        % Construct the full name for this field
        currentName = strcat(parentName, '.', fields{i});
        value = structure.(fields{i});
        
        % If the value is another structure, recurse into it
        if isstruct(value)
            if isequal(size(value), [1, 1])
                % Recurse into 1x1 struct only
                results_row = addNestedFieldsToResultsRow(value, currentName, results_row);
            end
        elseif isscalar(value)  % Only add if the field is a scalar (1x1)
            % For non-structure scalar fields, add them to results_row with a valid name
            valid_field_name = matlab.lang.makeValidName(currentName);
            results_row.(valid_field_name) = value;
        end
    end
end

