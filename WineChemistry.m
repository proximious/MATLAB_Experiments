function WineChemistry()
    [X, y, data_table] = linear_fitter_read_file();

    names = data_table.Properties.VariableNames(1:end-1);

    [selected_X, ~] = removeHighVIF(X, 5, names);

    [A, ~] = linear_fitter(selected_X, y);
    A_table = table(A(:), 'VariableNames', {'Coefficients A'});
    disp(A_table);
    
    R2_values = find_R2_values(selected_X);
    R2_table = table(R2_values(:), 'VariableNames', {'R Squared Values'});
    disp(R2_table);

end

function [X, y, data_table] = linear_fitter_read_file()
    data_table = readtable("shared\winequality-red.csv", "VariableNamingRule", "preserve");
    
    X = data_table(:, 1:end-1);  % All columns except the last one (x values)
    y = data_table(:, end);      % Last column is the output (y values)
    
    X = X{:,:};
    y = y{:,:};
end

function [A, X] = linear_fitter(X, y)
    % Add a column of ones as the first column for the bias term
    X = [ones(size(X, 1), 1), X]; 

    A = pinv(X' * X) * X' * y;
end

function vif_values = computeVIF(X)   
    [~, p] = size(X);  % n: number of samples, p: number of predictors
    vif_values = zeros(p, 1);  % Preallocate VIF array
    
    for i = 1:p
        Xi = X(:, i);                % Isolate the i-th predictor
        X_rest = X(:, [1:i-1, i+1:p]); % All other predictors
        
        % Fit linear regression model: Xi ~ X_rest
        b = X_rest \ Xi;    % Regression coefficients
        Xi_pred = X_rest * b;  % Predicted values
        R2 = 1 - sum((Xi - Xi_pred).^2) / sum((Xi - mean(Xi)).^2);  % R²
        
        vif_values(i) = 1 / (1 - R2);  % Compute VIF
    end
end

function R2_values = find_R2_values(X)
    [~, p] = size(X);  % n: number of samples, p: number of predictors
    R2_values = zeros(p, 1);  % Preallocate R2 array

    for i = 1:p
        Xi = X(:, i);                % Isolate the i-th predictor
        X_rest = X(:, [1:i-1, i+1:p]); % All other predictors
        
        % Fit linear regression model: Xi ~ X_rest
        b = X_rest \ Xi;    % Regression coefficients
        Xi_pred = X_rest * b;  % Predicted values
        R2 = 1 - sum((Xi - Xi_pred).^2) / sum((Xi - mean(Xi)).^2);  % R²
        
        R2_values(i) = R2;
    end
end

function [selected, removed] = removeHighVIF(X, threshold, names)

    removed = []; % Store removed variables
    iteration = 1;
    
    while true
        vif_values = computeVIF(X); % Compute VIF for current X
        fprintf('\nIteration %d - VIF Analysis:\n', iteration);
        
        % Display VIF scores
        for i = 1:length(vif_values)
            fprintf('Variable %s: VIF = %.2f\n', names{1, i}, vif_values(i));
        end
        
        % Find the variable with the highest VIF above the threshold
        [maxVIF, maxIdx] = max(vif_values);
        
        % Stop when all VIF values are below the threshold
        if maxVIF <= threshold
            break;
        end
        
        % Remove the variable with the highest VIF
        fprintf('Removing Variable %s (VIF = %.2f)\n', names{1, maxIdx}, maxVIF);
        
        % Store removed variable index BEFORE removing it, otherwise the
        % indexing is bad
        removed = [removed, names{1, maxIdx}]; 
        
        % remove the column with the removed value
        X(:, maxIdx) = [];
        
        % removing the associated name
        names(maxIdx) = [];
        iteration = iteration + 1;
    end 
        if size(removed) == 0
            disp("No removed values");
        else
            % display the removed values
            fprintf('Removed values: ');
            disp(removed);
        end
       
        % display the table
        T = table(names(:), vif_values(:), 'VariableNames', {'Variable Name', 'VIF Value'});
        disp(T);

        selected = X;
end
