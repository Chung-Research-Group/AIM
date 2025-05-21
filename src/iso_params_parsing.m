function opt = iso_params_parsing(isotherm_model, fitted_iso_params)
    % Function to parse the isotherm model parameters after fitting for
    % tabular representation
    switch(isotherm_model)
        case 'SS-Langmuir'
            opt.name_params = {'q_sat', 'b'};
            opt.unit_params = {'(mol/kg)', '(1/Pa)'};
            opt.express = {'$q^{*}(P) = \frac{q_{sat}bP}{1 + bP}$'};    
        case 'DS-Langmuir'
            opt.name_params = {'q_b_sat', 'b', 'q_d_sat', 'd'};
            opt.unit_params = {'(mol/kg)', '(1/Pa)', '(mol/kg)', '(1/Pa)'};
            opt.express = {'$q^{*}(P) = \frac{q_{b,sat}bP}{1 + bP} + \frac{q_{d,sat}dP}{1 + dP}$'};
        case 'SS-Langmuir-Freundlich'
            opt.name_params = {'q_sat', 'b', 'n'};
            opt.unit_params = {'(mol/kg)', '(1/Pa^n)', '(-)'};
            opt.express = {'$q^{*}(P) = \frac{q_{sat}bP^{n}}{1 + bP^{n}}$'};
        case 'DS-Langmuir-Freundlich'
            opt.name_params = {'q_b_sat', 'b', 'n_b', 'q_d_sat', 'd', 'n_d'};
            opt.unit_params = {'(mol/kg)', '(1/Pa^n_b)', '(-)', '(mol/kg)', '(1/Pa^n_d)', '(-)'};
            opt.express = {'$q^{*}(P) = \frac{q_{b,sat}bP^{n_{b}}}{1 + bP^{n_{b}}} + \frac{q_{d,sat}dP^{n_{d}}}{1 + dP^{n_{d}}}$'};
        case 'Quadratic'
            opt.name_params = {'q_sat', 'b', 'c'};
            opt.unit_params = {'(mol/kg)', '(1/Pa)', '(1/Pa^2)'};
            opt.express = {'$q^{*}(P) = q_{sat}\left(\frac{bP+2cP^{2}}{1 + bP + cP^{2}}\right)'};
        case 'Temkin'
            opt.name_params = {'q_sat', 'b', 'theta'};
            opt.unit_params = {'(mol/kg)', '(1/Pa)', '(-)'};
            opt.express = {'$q^{*}(P) = q_{sat}\left(\frac{bP}{1 + bP}\right) \\+ q_{sat}\theta \left(\frac{bP}{1 + bP}\right)^{2} \left(\frac{bP}{1 + bP}- 1\right)$'};
        case 'BET'
            opt.name_params = {'q_sat', 'b', 'c'};
            opt.unit_params = {'(mol/kg)', '(1/Pa)', '(1/Pa)'};
            opt.express = {'$q^{*}(P) = \frac{q_{sat}bP}{(1 -cP)(1 - cP + bP)}$'};
        case 'Sips'
            opt.name_params = {'q_sat', 'b', 'n'};
            opt.unit_params = {'(mol/kg)', '(1/Pa)', '(-)'};
            opt.express = {'$q^{*}(P) = \frac{q_{sat}(bP)^{\frac{1}{n}}}{1 + (bP)^{\frac{1}{n}}}$'};
        case 'Toth'
            opt.name_params = {'q_sat', 'b', 'n'};
            opt.unit_params = {'(mol/kg)', '(1/Pa)', '(-)'};
            opt.express = {'$q^{*}(P) = \frac{q_{sat}bP}{(1 + (bP)^{n})^{\frac{1}{n}}}$'};
        case 'Dubinin-Astakhov'
            opt.name_params = {'q_sat', 'K', 'n'};
            opt.express = '$$x=\frac{P}{P_{0}},\quad q^{*}(x) = q_{sat}e^{-\left(\frac{1}{K}\times ln\frac{1}{x}\right)^{n}}$$';
        case 'Klotz'
            opt.name_params = {'q_sat', 'K', 'C', 'n'};
            opt.express = '$$x=\frac{P}{P_{0}},\quad s=Kx\\ q^{*}(s) = q_{sat}\frac{Cs\{1-(1+n)s^{n}+ns^{n+1}\}}{(1-s)\{1+(C-1)s-Cs^{n+1}\}}$$';
        case 'Do-Do'
            opt.name_params = {'q_sat', 'f', 'K_1', 'K_2', 'alpha', 'beta'};
            opt.express = ['$$x=\frac{P}{P_{0}}\\ q^{*}(x) = q_{sat}\left(f\frac{K_{1}x\{1-(1+\beta)x^{\beta}+\beta x^{\beta +1}\}}{(1-x)\{1+(K_{1}-1)x-K_{1}x^{\beta +1}\}}'...
                                        +'\\+(1-f)\frac{K_{2}x^{\alpha}}{1+K_{2}x^{\alpha}$$\right)'];
        case 'Structural-Transition-Adsorption'
            opt.name_params = {'q_NP_sat', 'b_NP', 'q_LP_sat', 'b_LP', 's', 'P_tr'};
            opt.express = ['$$y(P)=\left(\frac{1+b_{NP}P_{tr}}{1+b_{NP}P}\right)^{q_{NP,sat}}'...
                            '\times\left(\frac{1+b_{LP}P}{1+b_{LP}P_{tr}}\right)^{q_{LP,sat}},'...
                            '\\\sigma(P)=\frac{y^{s}}{1+y^{s}},' ...
                            '\\q^{*}(P) = \left(1-\sigma\right)\left(\frac{q_{NP,sat}b_{NP}P}{1+b_{NP}P}\right)'...
                            '+\sigma\left(\frac{q_{LP,sat}b_{LP}P}{1+b_{LP}P}\right)'];
        
        case 'Virial'
            opt.name_params = {};
            opt.express = ['$$\ln P = \ln q + \frac{1}{T}\sum_{i=0}^{m} a_{i}\;q^{i}'...
                            '+\sum_{j=0}^{n} b_{j}\;q^{j}'];

        otherwise
            error('Unrecognized isotherm model %s', isotherm_model);
    end
    
    if ~(isempty(fitted_iso_params)) && ~strcmpi(isotherm_model, 'Virial')
        if ~(numel(opt.name_params) == length(fitted_iso_params))
            error("Mismatch of number of parameters and their values...!")
        end
        opt.num_params = length(fitted_iso_params);
    else
        opt.num_params = length(opt.name_params);
    end
end