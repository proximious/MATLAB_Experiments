function linear_fitter_test()
    manufacture_data()
end

function A = linear_fitter(data)
    % Extract X (input variables) and y (output variable)
    X = data(:, 1:end-1);  % All columns except the last one (x values)
    y = data(:, end);      % Last column is the output (y values)

    % Add a column of ones as the first column for the bias term
    X = [ones(size(X, 1), 1), X]; 

    A = pinv(X' * X) * X' * y;
end

function manufacture_data()
    % Step 1: Generate fake data
    number_of_points = 100;
    x = linspace(-10, 10, number_of_points)';  % Generate 100 x values from -10 to 10

    % found an example online
    % y = 4 - 0.5x + 0.05x^2 + 1.2cos(x) + noise
    true_params = [4, -0.5, 0.05, 1.2];  % [a0, a1, a2, a3]
    noise = 1 + 2 * randn(number_of_points, 1);  % Random noise with mean 1, std dev 2
    y = true_params(1) + (true_params(2) * x) + (true_params(3) * x.^2) + (true_params(4) * cos(x)) + noise;

    % Step 2: Construct Data Matrix
    X = [x, x.^2, cos(x)];  % Create feature matrix
    data = [X, y];  % Append output y as last column

    % Step 3: Compute Best-Fit Parameters using linear_fitter
    fitted = linear_fitter(data);
    
    % Step 4: Plot Original Data vs. Fitted Curve
    y_fitted = fitted(1) + (fitted(2) * x) + (fitted(3) * x.^2) + (fitted(4) * cos(x));

    figure;
    scatter(x, y, 'blue', 'filled'); 
    hold on;  % Scatter plot of noisy data
    plot(x, y_fitted, 'red', 'LineWidth', 2);  % Plot fitted curve
    title('Linear Fit vs. Noisy Data');
    legend('Noisy Data', 'Fitted Curve');
    grid on;
    
    % Step 5: Display Results
    disp('True Parameters:');
    disp(true_params(:));  % Display true parameter vector
    disp('Fitted Parameters:');
    disp(fitted);  % Display estimated parameters
end