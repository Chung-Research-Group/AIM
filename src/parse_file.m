function file_data = parse_file(filename)
    % Function to parse *.csv or *.dat or *.txt files containing isotherm
    % data
    
    % Initialize the struct object
    file_data = struct();

    file_data.meta_data = dictionary;
    file_data.data = [];
    
    % Matching Pattern
    match_pat = regexpPattern('\-?(\d+(\.\d*)?|\.\d+)');
    
    % Opening the file
    fileID = fopen(filename, "r");
    
    % Number of meta data entries
    num_meta = 0;
    num_data = 0;

    % Place holder for isotherm data
    isotherm_data = NaN(1000, 2);

    %% Looping every line of file
    while ~feof(fileID)        
        current_line = fgetl(fileID);
        
        % Ignoring empty line
        if isempty(current_line)
           continue
        end

        % Check for char data type
        if ischar(current_line)
            
            % Trimming whitespaces
            current_line = strip(current_line);
            
            % Removing any quotation marks or extra commas
            current_line = replace(current_line, {'"', ''''}, '');

            %% Handling meta data
            if startsWith(current_line, '#')
                handle_meta_data();
                continue
            end
    
            %% Handling numerical data
            if startsWith(current_line, match_pat)
                handle_numeric_data();
            end
        end
    end

    % Closing the file
    fclose(fileID);
    
    % isotherm_data = readmatrix(filename);
    % file_data.data = isotherm_data;

    file_data.data = isotherm_data(1:num_data, :);
    %
    %% Miscellenous Functions
    % Function to handle meta data
    function [] = handle_meta_data()
        num_meta = num_meta + 1;
        try
            current_line = replace(current_line, ',', '');
            tokens = split(string(current_line));
            key = tokens(1);
            value = join(tokens(2:end), ' ');
            file_data.meta_data(key) = value;
        catch
            return
        end
    end
    %
    
    % Function to handle numeric data
    function [] = handle_numeric_data()
        try
            num_data = num_data + 1;
            current_line = strip(current_line, "both", ",");
            current_line = string(current_line);
            current_line = replace(current_line, ...
                            [",", "|", "\t", "\n", "`", "'", "\",...
                            "/", newline], " ");
            tokens = split(current_line);
            tokens(tokens == "") = [];
            data = reshape(str2double(tokens), 1, []); 
            isotherm_data(num_data, 1:2) = data(1, 1:2);
            
        catch
            isotherm_data(num_data, 1:2) = [NaN, NaN];
        end
    end
end