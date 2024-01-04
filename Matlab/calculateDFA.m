function dfa = calculateDFA(ts)
    v_gluc = true(1, length(ts)); % a "shadow" time series is created to account for NAs
    missing_val = find(isnan(ts));
    v_gluc(missing_val) = false; % identifying missing values
    v_types = repmat('B', 1, length(ts));  % identifying B (present), A (missing), PP (pre-NA segment), PF (post-NA segment)
    v_types(missing_val) = 'A';

    PP = [];
    PF = [];
    for i = missing_val
        if ~isnan(ts(i - 1))
            PP = [PP, i - 1];
        end
        if ~isnan(ts(i + 1))
            PF = [PF, i + 1];
        end
    end

    if ~isempty(PP)
        for j = 1:length(PP)
            fragm = linspace(ts(PP(j)), ts(PF(j)), PF(j) - PP(j) + 1);
            ts(PP(j):PF(j)) = round(fragm); % segment interpolation
        end
    end

    % Integration
    integrat = false;
    if integrat
        ts = ts - mean(ts);
        for int = 2:length(ts)
            ts(int) = ts(int) + ts(int - 1);
        end
    end

    w_size = [3, 4, 6, 8, 9, 12, 16, 18, 24, 32, 36, 48, 72, 96, 144, 288]; % window sizes
    tol = 0.2; % "tolerance" to interpolation
    fn_tot = [];

    for n1 = w_size
        dkblock = 0;
        notol = 0; % counter for windows >tol
        init = 1:n1:max(w_size); % start point
        n_vent = length(init);

        for n2 = init
            block = n2:(n2 + n1 - 1); % segment on analysis
            if sum(v_gluc(block)) / length(block) > tol % tol condition == T; detrending begins
                sery = ts(block); % segment selection
                str_line = fitlm(block', sery); % linear regression
                coeff = str_line.Coefficients.Estimate(2);
                inter = str_line.Coefficients.Estimate(1);
                y_teor = inter + (coeff * block');
                dk = sum((sery - y_teor').^2); % segment detrending
                dkblock = dkblock + dk;
            else % tol condition == F
                notol = notol + 1;
            end
        end

        meandk = dkblock / (n_vent - notol); % average of detrended segments
        dkblock = meandk * n_vent; % result including notol windows
        denom_fn = n_vent * n1;
        fn_bl = sqrt((1 / denom_fn) * dkblock);
        fn_tot = [fn_tot, fn_bl]; % Fn for the entire time series
    end

    calc_regress = fitlm(log(w_size)', log(fn_tot));
    plot(log(w_size), log(fn_tot), 'o', 'MarkerSize', 8);
    xlabel('log(window)');
    ylabel('log(Fn)');
    hold on;
    plot(log(w_size), calc_regress.Fitted, 'LineWidth', 2);
    legend('Data Points', 'Linear Regression');
    hold off;

    dfa = calc_regress.Coefficients.Estimate(2);
    result = sprintf('DFA = %.4f', dfa);
    disp(result);


end
