function aif_data = parse_aif(filename)
    % Function to parse AIF files
    
    % Initialize the struct object
    aif_data = struct();
    
    aif_data.comments = {};
    aif_data.meta_data = dictionary;
    aif_data.loop_data = [];
    aif_data.loop_keys = {};
    
    % Variables to identify loops in aif file
    loop_flag = 0;
    loop_count = 0;
    
    % Matching Expression
    match_exp = '\-?(\d+(\.\d*)?|\.\d+)';

    % Opening the file
    fileID = fopen(filename, "r");

    %% Looping every line of file
    while ~feof(fileID)
        
        current_line = fgetl(fileID);
    
        % Trimming whitespaces
        current_line = strip(current_line);
    
        % Ignoring empty line
        if isempty(current_line)
           continue
        end
    
        %% Handling Comments, comments will be stored as cell array with strings
        if startsWith(current_line, '#')
            aif_data.comments{end+1} = {string(current_line)};
            continue
        end
        
        %% Checking for loop keyword
        if startsWith(current_line, 'loop_')
            loop_count = loop_count + 1;

            % Only one loop per aif will be read
            if loop_count > 1
                loop_flag = 0;
            else
                loop_flag = 1;
                loop_values = [];
                continue
            end
        end
    
        %% Parsing meta-data        
        if startsWith(current_line, '_') && ~loop_flag
            tokens = split(string(current_line));
            key = tokens(1);
            value = join(tokens(2:end), ' ');
            aif_data.meta_data(key) = value;
            continue
        end
        
        %% Parsing loop data
        if loop_flag
            % Parsing column names specified for loop columns
            if startsWith(current_line, '_')
                aif_data.loop_keys(end+1) = {char(current_line)};
                continue

            % Parsing loop numeric data
            else
                values = regexp(current_line, match_exp, "match");           
                values = double(string(values));
                if all(~isnan(values))
                    try
                        loop_values = vertcat(loop_values, values);
                    catch ME1
                        error("Failed to concatenate loop data in AIF file, %s", ME1.message);
                    end
                end 
                continue
            end
        end
    end

    % Closing the file
    fclose(fileID);

    aif_data.loop_data = loop_values;
end