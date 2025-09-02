function total_fun()
    % given values
    c = 0.05;
    L = 20;
    N = 10;
    K = 100; % number of divisions

    % max and tolerance to find out when to end looping
    tol = 1e-6; % convergence tolerance
    max_iter = 100;

    lambda = 1 / (1 + L) - 0.001;

    % making a constant to use 
    c1 = sqrt(1/(3*c));

    % defining our functions given from the writeup
    B = @(lambda) c1 * ...
        trapezoider(@(x) sqrt(max(0, 1./(1 + x) - lambda)), 0, L, K) - N;

    B_prime = @(lambda) -0.5 * c1 * ...
        trapezoider(@(x) 1 ./ sqrt(max(0, 1./(1 + x) - lambda)), 0, L, K);

    % performing newton-raphson method
    % in this case, 
    %   x_n = lambda
    %   x_n+1 = lambda_new
    %   f(x_n) = b_val
    %   f'(x_n) = b_prime_val
    for iter = 1:max_iter
        b_val = B(lambda);
        b_prime_val = B_prime(lambda);
        
        % Update
        lambda_new = lambda - (b_val / b_prime_val);
        
        % Check convergence
        if abs(lambda_new - lambda) < tol
            lambda = lambda_new;
            break;
        end
    
        lambda = lambda_new;
    end
    
    fprintf('Converged lambda: %.6f (in %d iterations)\n', lambda, iter);
    fprintf('Expected rho = 0 starting around x = %.2f\n', 1/lambda - 1);

    x_vals = linspace(0, L, K+1);
    rho_vals = sqrt(1/(3*c) * (1./(1 + x_vals) - lambda));
    
    % Set complex values to 0
    rho_vals(~isreal(rho_vals)) = 0;
    
    % Plot
    figure;
    plot(x_vals, rho_vals, 'b-', 'LineWidth', 2);
    xlabel('x');
    ylabel('\rho(x)');
    title('Optimal Crowd Density \rho(x)');
    grid on;
end


function [A] = trapezoider(f, x_L, x_U, K)
    % following https://byjus.com/maths/trapezoidal-rule/

    % first have to get the step size
    % delta x = (b - a) / N
    delta_x = (x_U - x_L) / K;

    % get all the x values that we need
    x_values = linspace(x_L, x_U, K+1);

    % get all the y values that we need
    y_values = f(x_values);

    % Set values that are not real to zero
    y_values(~isreal(y_values)) = 0;

    % apply the trapezoidal rule
    A = (delta_x/2) * (y_values(1) + 2*sum(y_values(2:end-1)) + y_values(end));
end