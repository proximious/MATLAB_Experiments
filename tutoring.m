%%%
% run examples
%%%
function tutoring()
    convergenceExample()
    fprintf("------\n")
    roundOffExample()
    fprintf("------\n")
    truncationExample()
    fprintf("------\n")
    
    OneNewtonRaphsonExample()
    fprintf("------\n")
    MultiNewtonRaphsonExample()
    fprintf("------\n")
    
    LUDecompositionExample()
    fprintf("------\n")

    secantMethodExample()
    fprintf("------\n")        
end

%%%
% Convergence, Round-off, Truncation
%%%
function convergenceExample()
    fprintf("Convergence Example\n")

    f = @(x) x^2 - 4;       % Function: f(x) = x^2 - 2
    df = @(x) 2*x;          % Derivative: f'(x) = 2x
    
    x = 1;                  % Initial guess
    tol = 1e-6;             % Tolerance level
    max_iter = 50;          % Maximum iterations
    iter = 0;
    
    while abs(f(x)) > tol && iter < max_iter
        x = x - f(x) / df(x); % Newton-Raphson update
        iter = iter + 1;
        fprintf('Iteration %d: x = %.6f, Error = %.6e\n', iter, x, abs(f(x)));
    end
    
    disp(['Root approximation: ', num2str(x)]);
end

function roundOffExample()
    fprintf("Round-off Example\n")

    x_true = 1/3;            % Exact value
    x_approx = round(x_true, 4); % Approximate to 4 decimal places
    
    abs_error = abs(x_true - x_approx); % Absolute error
    rel_error = abs_error / abs(x_true); % Relative error
    
    fprintf('True value: %.10f\n', x_true);
    fprintf('Approx value: %.4f\n', x_approx);
    fprintf('Absolute error: %.4e\n', abs_error);
    fprintf('Relative error: %.4e\n', rel_error);
end

function truncationExample()
    fprintf("Truncation Example\n")

    x = pi / 4; % True value (45 degrees)
    
    % True sine value
    true_sin = sin(x);
    
    % Taylor series approximation (3 terms)
    approx_sin = x - (x^3) / factorial(3) + (x^5) / factorial(5);
    
    % Calculate truncation error
    truncation_error = abs(true_sin - approx_sin);
    
    fprintf('True sin(pi/4) = %.10f\n', true_sin);
    fprintf('Approx sin(pi/4) = %.10f\n', approx_sin);
    fprintf('Truncation error = %.10e\n', truncation_error);
end

%%%
% Newton-Raphson Methods
%%%

function OneNewtonRaphsonExample()
    % Define function and derivative
    f = @(x) cos(x);  
    df = @(x) -sin(x);    
    
    % Starting value
    x = 4;  
    tol = 1e-6;  
    max_iter = 50;  
    iter = 0;
    
    % Newton-Raphson iteration
    while abs(f(x)) > tol && iter < max_iter
        x = x - f(x) / df(x);  % Newton-Raphson formula
        iter = iter + 1;
        fprintf('Iteration %d: x = %.6f, Error = %.6e\n', iter, x, abs(f(x)));
    end
    
    disp(['Root approximation: ', num2str(x)]);
end

function MultiNewtonRaphsonExample()
    % Define the nonlinear system
    F = @(x) [x(1)^2 + x(2)^2 - 4;
            x(1) - sin(x(2))];

    % Define the Jacobian matrix
    J = @(x) [2*x(1), 2*x(2);
               1, -cos(x(2))];
    
    % Initial guess
    x = [1; 1];  % Column vector
    tol = 1e-6;
    max_iter = 50;
    iter = 0;
    
    while norm(F(x)) > tol && iter < max_iter
        dx = -J(x) \ F(x);  % Solve J(x) * dx = -F(x)
        x = x + dx;
        iter = iter + 1;
        fprintf('Iteration %d: x = %.6f, y = %.6f, Error = %.6e\n', iter, x(1), x(2), norm(F(x)));
    end
    
    disp(['Solution: x = ', num2str(x(1)), ', y = ', num2str(x(2))]);

end

%%%
% LU Decomposition
%%%

function LUDecompositionExample()
    % Define matrix A and vector b
    A = [2 -1  1;
         3  3  9;
         3  3  5];
    
    b = [1; 6; 4];
    
    % Perform LU decomposition
    [L, U] = lu(A);  % MATLAB’s built-in LU decomposition
    
    % Solve Ld = b using forward substitution
    d = L \ b;
    
    % Solve Ux = d using backward substitution
    x = U \ d;
    
    % Display results
    disp('LU decomposition:');
    disp('L ='), disp(L);
    disp('U ='), disp(U);
    disp('Solution x ='), disp(x);
end


%%%
% Secant method
%%%
function secantMethodExample()
    f = @(x) x^2 - 4;      
    % Starting value
    x0 = 1;
    x1 = 3;  
    tol = 1e-6;  
    max_iter = 50;  
    iter = 0;
    error = abs(x1 - x0);

    while error > tol && iter < max_iter
        iter = iter + 1;  % Increment iteration counter
        
        % Compute function values
        f_x0 = f(x0);
        f_x1 = f(x1);
        
        % Check for zero denominator (avoid division by zero)
        if abs(f_x1 - f_x0) < eps
            error('Denominator close to zero. Choose different initial guesses.');
        end
        
        % Compute next approximation using the Secant formula
        x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);
        
        % Update error estimate
        error = abs(x2 - x1);

        % Display iteration results
        fprintf('Iteration %d: x = %.6f, y = %.6f, Error = %.6e\n', iter, x2, f(x2), error);
        
        % Check stopping condition
        if error < tol || abs(f(x2)) < tol
            fprintf("------\n")
            fprintf('Solution: x =  %.6f\n', x2);
            return;
        end
        
        % Update for next iteration
        x0 = x1;
        x1 = x2;
    end
    
    % If max iterations reached without convergence
    fprintf("------\n")
    error('Secant method did not converge within the maximum iterations.');
end



