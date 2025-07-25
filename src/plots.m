function [] = plots(ax , L, N, variable_to_plot, solution, plot_type, time_value, parameter_set)    
    %% Prameters unpacking  
    comp_names = parameter_set.CompNames;
    comp_num = parameter_set.CompNum;
    params = ProcessInputParameters(parameter_set);
    inlet_MF = [params(24:27); 1-sum(params(24:27))]; 
    inlet_MF(inlet_MF==0) = 1e-05;

    dz = L/N ;
    t = solution.t(:, 1) ;
    nodes = (1:1:N+2)';
    len_column = zeros(size(nodes));
    len_column(2:N+1) = nodes(1:N) .* dz - dz/2;
    len_column(1) = len_column(2) - dz/2;
    len_column(end) = len_column(N+1) + dz/2;
%     time_value_index = round(time_value/0.01 + 1);
    time_value_index = round(time_value/1 + 1);
    
    line_width = 2.5;
    scatter_marker_size = 3.0;

    fontsize = 20;
    fontname = 'Calibri';
    colororder(ax, "gem12");
    % ax.FontSize = 18;
    %
    %% Pressure 
    if strcmpi(variable_to_plot, "Pressure")
        P = solution.P ./ 1e3;  % Conversion to kPa for better visualization
        % Across the column plot
        if strcmpi(plot_type, "Across the column")
            p = plot(ax, len_column, P(time_value_index,:), '-ko');
            p.LineWidth = line_width;
%             p.MarkerSize = 2.0;

            

            ax.Title.String = "Pressure Profile across the column";
            xlabel(ax, 'Length (m)', FontSize = fontsize, FontWeight="bold", FontName='AvantGrande');
        
        % At the Outlet plot    
        else
            p = plot(ax, t, P(:, end), '-ko');
            p.LineWidth = line_width;
            p.MarkerSize = scatter_marker_size;
            p.MarkerIndices = 1:1500:length(t);

            ax.Title.String = "Pressure Profile at column outlet";
            xlabel(ax, 'Time (s)', FontSize = fontsize, FontWeight="bold", FontName='AvantGrande');
        end

        ax.Box = "off";
        ax.FontName = "Arial";
        ax.Color = [1 1 1];
        
        ax.Title.FontWeight ="bold";
        ax.Title.FontSize = fontsize;
        ax.Title.FontName = fontname;
        ylabel(ax, 'Pressure (kPa)', FontSize = fontsize, FontWeight="bold", FontName='AvantGrande');
    end
%
    %% Temperature 
    if strcmpi(variable_to_plot, "Temperature")
        T = solution.T;
        % Across the column plot
        if strcmpi(plot_type, "Across the column")
            p = plot(ax, len_column, T(time_value_index,:), '-o');
            p.LineWidth = line_width;
%             p.MarkerSize = 2.0;


            ax.Title.String = "Temperature Profile across the column";
            xlabel(ax, 'Length (m)', FontSize = fontsize, FontWeight="bold", FontName='AvantGrande');
        
        % At the Outlet plot    
        else
            p = plot(ax, t, T(:, end), '-o');
            p.LineWidth = line_width;
            p.MarkerSize = scatter_marker_size;
            p.MarkerIndices = 1:1500:length(t);

            ax.Title.String = "Temperature Profile at column outlet";
            xlabel(ax, 'Time (s)', FontSize = fontsize, FontWeight="bold", FontName='AvantGrande');
            
        end

        ax.Box = "off";
        ax.FontName = "Arial";
        ax.Color = [1 1 1];
        p.LineWidth = line_width;
        p.Color = [0.85, 0.325, 0.098];
        ax.Title.FontWeight ="bold";
        ax.Title.FontSize = fontsize;
        ax.Title.FontName=fontname;
        % c = "Temperature ("+char(176)+"C)";
        c = "Temperature (K)";
        ylabel(ax, c, FontSize = fontsize, FontWeight="bold", FontName='AvantGrande');
    end
%
%%  Wall Temperature 
    if strcmpi(variable_to_plot, "Wall Temperature")
        Tw = solution.Tw;
        % Across the column plot
        if strcmpi(plot_type, "Across the column")
            p = plot(ax, len_column, Tw(time_value_index,:), '-o');
            p.LineWidth = line_width;

            ax.Title.String = "Wall Temperature Profile across the column";
            xlabel(ax, 'Length (m)', FontSize = fontsize, FontWeight="bold", FontName='AvantGrande');
        
        % At the Outlet plot    
        else
            p = plot(ax, t, Tw(:, end), '-o');
            p.LineWidth = line_width;
            p.MarkerSize = scatter_marker_size;
            p.MarkerIndices = 1:1500:length(t);

            ax.Title.String = "Wall Temperature Profile at column outlet";
            xlabel(ax, 'Time (s)', FontSize = fontsize, FontWeight="bold", FontName='AvantGrande');
        end

        ax.Box = "off";
        ax.FontName = "Arial";
        ax.Color = [1 1 1];
        p.LineWidth = line_width;
        p.Color = [0.65, 0.625, 0.198];
        ax.Title.FontWeight ="bold";
        ax.Title.FontSize = fontsize;
        ax.Title.FontName=fontname;
        % c = "Temperature ("+char(176)+"C)";
        c = "Wall Temperature (K)";
        ylabel(ax, c, FontSize = fontsize, FontWeight="bold", FontName='AvantGrande');
    end
%
    %% Breakthrough Profiles
    if strcmpi(variable_to_plot, "Mole Fractions")
        C1 = solution.C1;
        C2 = solution.C2;
        C3 = solution.C3;
        C4 = solution.C4;
        C5 = solution.C5;
    

        if strcmpi(plot_type, "Across the column")
            C1_col = C1(time_value_index, :);
            C2_col = C2(time_value_index, :);
            C3_col = C3(time_value_index, :);
            C4_col = C4(time_value_index, :);
            C5_col = C5(time_value_index, :);
    
            % Plotting
            if comp_num == 2
                p = plot(ax, len_column, C5_col, '-go', len_column, C1_col, '-mo');
            elseif comp_num==3
                p = plot(ax, len_column, C5_col, '-go', len_column, C1_col, '-mo', len_column, C2_col, "-o");
            elseif comp_num==4
                p = plot(ax, len_column, C5_col, '-go', len_column, C1_col, '-mo', len_column, C2_col, "-o", len_column, C3_col, "-bo");
            elseif comp_num==5
                p = plot(ax, len_column, C5_col, '-go', len_column, C1_col, '-mo', len_column, C2_col, "-o", len_column, C3_col, "-bo", len_column, C4_col, "-ro");
            end

            title(ax, "Composition Profiles across the column")
            xlabel(ax, "Length (m)", FontSize = fontsize, FontName='AvantGrande', FontWeight="bold");
            
        % Outlet Dry based MF calculation
        else
            C1_out = C1(:, end);
            C2_out = C2(:, end);
            C3_out = C3(:, end);
            C4_out = C4(:, end);
            C5_out = C5(:, end);
    
            % Plotting
            if comp_num == 2
                p = plot(ax, t, C5_out, '-go', t, C1_out, '-mo');
            elseif comp_num == 3
                p = plot(ax, t, C5_out, '-go', t, C1_out, '-mo', t, C2_out, "-o");
            elseif comp_num == 4
                p = plot(ax, t, C5_out, '-go', t, C1_out, '-mo', t, C2_out, "-o", t, C3_out, "-bo");
            elseif comp_num == 5
                p = plot(ax, t, C5_out, '-go', t, C1_out, '-mo', t, C2_out, "-o", t, C3_out, "-bo", t, C4_out, '-ro');            
            end

            for i = 1:comp_num
                p(i).MarkerSize = scatter_marker_size;
                p(i).MarkerIndices = 1:1500:length(t);
            end
            title(ax, "Composition Profiles at the outlet")
            xlabel(ax, "Time (s)", FontSize = fontsize, FontName='AvantGrande', FontWeight="bold");
        end
        

        for i = 1:comp_num
            p(i).LineWidth = line_width;
        end
        
        ax.Box = "off";
        ax.FontName = "Arial";
        ax.Color = [1 1 1];
        ax.YLim = [0 1];

        ax.YTick = [0 0.2 0.4 0.6 0.8 1.0];
        ax.Title.FontWeight ="bold";
        ax.Title.FontSize = fontsize;
        ax.Title.FontName=fontname;
        % names = [comp_names(comp_num), comp_names(1:comp_num-1)];
        names = comp_names;
        legend(ax, convertStringsToChars(names))
        ylabel(ax, "Mole fraction", FontSize = fontsize, FontName='AvantGrande', FontWeight="bold");
    end
%
    %% Molar Loading
    if strcmpi(variable_to_plot, "Molar Loadings")
        x_C1 = solution.x1C1;
        x_C2 = solution.x2C2;
        x_C3 = solution.x3C3;
        x_C4 = solution.x4C4;
        x_C5 = solution.x5C5;

        if strcmpi(plot_type, "Across the column")
            C1_col = x_C1(time_value_index, :);
            C2_col = x_C2(time_value_index, :);
            C3_col = x_C3(time_value_index, :);
            C4_col = x_C4(time_value_index, :);
            C5_col = x_C5(time_value_index, :);
    
            % Plotting
            if comp_num == 2
                p = plot(ax, len_column, C5_col, '-go', len_column, C1_col, '-mo');
            elseif comp_num==3
                p = plot(ax, len_column, C5_col, '-go', len_column, C1_col, '-mo', len_column, C2_col, "-o");
            elseif comp_num==4
                p = plot(ax, len_column, C5_col, '-go', len_column, C1_col, '-mo', len_column, C2_col, "-o", len_column, C3_col, "-bo");
            elseif comp_num==5
                p = plot(ax, len_column, C5_col, '-go', len_column, C1_col, '-mo', len_column, C2_col, "-o", len_column, C3_col, "-bo", len_column, C4_col, "-ro");
            end

            title(ax, "Molar Loadings across the column")
            xlabel(ax, "Length (m)", FontSize = fontsize, FontName='AvantGrande', FontWeight="bold");
            
        % Outlet Dry based MF calculation
        else
            C1_out = x_C1(:, end);
            C2_out = x_C2(:, end);
            C3_out = x_C3(:, end);
            C4_out = x_C4(:, end);
            C5_out = x_C5(:, end);
    
            % Plotting
            if comp_num == 2
                p = plot(ax, t, C5_out, '-go', t, C1_out, '-mo');
            elseif comp_num == 3
                p = plot(ax, t, C5_out, '-go', t, C1_out, '-mo', t, C2_out, "-o");
            elseif comp_num == 4
                p = plot(ax, t, C5_out, '-go', t, C1_out, '-mo', t, C2_out, "-o", t, C3_out, "-bo");
            elseif comp_num == 5
                p = plot(ax, t, C5_out, '-go', t, C1_out, '-mo', t, C2_out, "-o", t, C3_out, "-bo", t, C4_out, '-ro');            
            end

            for i = 1:comp_num
                p(i).MarkerSize = scatter_marker_size;
                p(i).MarkerIndices = 1:1500:length(t);
            end
            title(ax, "Molar Loadings at the outlet")
            xlabel(ax, "Time (s)", FontSize = fontsize, FontName='AvantGrande', FontWeight="bold");
        end
    
        for i = 1:comp_num
            p(i).LineWidth = line_width;
        end
        ax.Box = "off";
        ax.FontName = fontname;
        ax.Color = [1 1 1];
        ax.Title.FontWeight ="bold";
        ax.Title.FontSize = fontsize;
        ax.Title.FontName= fontname;
        % names = [comp_names(comp_num), comp_names(1:comp_num-1)];
        names = comp_names;
        legend(ax, convertStringsToChars(names))
        ylabel(ax, "Molar Loadings (mol/kg)", FontSize = fontsize, FontName='AvantGrande', FontWeight="bold");
    end
    ax.FontName = fontname;
    ax.FontWeight = 'bold';

    %% Normalized MF
 
    if strcmpi(variable_to_plot, "Normalized Mole Fractions")
        C1 = solution.C1./inlet_MF(1);
        C2 = solution.C2./inlet_MF(2);
        C3 = solution.C3./inlet_MF(3);
        C4 = solution.C4./inlet_MF(4);
        C5 = solution.C5./inlet_MF(5);
    

        if strcmpi(plot_type, "Across the column")
            C1_col = C1(time_value_index, :);
            C2_col = C2(time_value_index, :);
            C3_col = C3(time_value_index, :);
            C4_col = C4(time_value_index, :);
            C5_col = C5(time_value_index, :);
    
            % Plotting
            if comp_num == 2
                p = plot(ax, len_column, C5_col, '-go', len_column, C1_col, '-mo');
            elseif comp_num==3
                p = plot(ax, len_column, C5_col, '-go', len_column, C1_col, '-mo', len_column, C2_col, "-o");
            elseif comp_num==4
                p = plot(ax, len_column, C5_col, '-go', len_column, C1_col, '-mo', len_column, C2_col, "-o", len_column, C3_col, "-bo");
            elseif comp_num==5
                p = plot(ax, len_column, C5_col, '-go', len_column, C1_col, '-mo', len_column, C2_col, "-o", len_column, C3_col, "-bo", len_column, C4_col, "-ro");
            end

            title(ax, "Composition Profiles across the column")
            xlabel(ax, "Length (m)", FontSize = fontsize, FontName='AvantGrande', FontWeight="bold");
            
        % Outlet Dry based MF calculation
        else
            C1_out = C1(:, end);
            C2_out = C2(:, end);
            C3_out = C3(:, end);
            C4_out = C4(:, end);
            C5_out = C5(:, end);
    
            % Plotting
            if comp_num == 2
                p = plot(ax, t, C5_out, '-go', t, C1_out, '-mo');
            elseif comp_num == 3
                p = plot(ax, t, C5_out, '-go', t, C1_out, '-mo', t, C2_out, "-o");
            elseif comp_num == 4
                p = plot(ax, t, C5_out, '-go', t, C1_out, '-mo', t, C2_out, "-o", t, C3_out, "-bo");
            elseif comp_num == 5
                p = plot(ax, t, C5_out, '-go', t, C1_out, '-mo', t, C2_out, "-o", t, C3_out, "-bo", t, C4_out, '-ro');            
            end

            for i = 1:comp_num
                p(i).MarkerSize = scatter_marker_size;
                p(i).MarkerIndices = 1:1500:length(t);
            end
            title(ax, "Composition Profiles at the outlet")
            xlabel(ax, "Time (s)", FontSize = fontsize, FontName='AvantGrande', FontWeight="bold");
        end
        

        for i = 1:comp_num
            p(i).LineWidth = line_width;
        end
        
        ax.Box = "off";
        ax.FontName = fontname;
        ax.Color = [1 1 1];

        ax.Title.FontWeight ="bold";
        ax.Title.FontSize = fontsize;
        ax.Title.FontName=fontname;
        % names = [comp_names(comp_num), comp_names(1:comp_num-1)];
        names = comp_names;
        legend(ax, convertStringsToChars(names))
        ylabel(ax, "Normalized Concentration y_i /y_{io}", FontSize = fontsize, FontName='AvantGrande', FontWeight="bold"); 
    end
    ax.FontName = fontname;
    ax.FontWeight = 'bold';
    ax.FontSize = fontsize;
end