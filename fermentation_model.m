function fermentation_model()
    % given parameters 
    a = 0.1; 
    b = 0.3; 
    c = 1.1; 
    d = 0.3; 
    g = 0.1; 
    x2_crit = 0.2;
    
    % Initial conditions
    guess = [0.8; 0.8; 0.2];
    t0 = 0;
    tf = 15;
    dt = 0.1;
    tspan = t0:dt:tf;

    seeds = [0.8; 0.8; 0.2];  % [s0, q0, y0] given seed values
    d_values = [0.01, 0.01, 0.01];
    
    target_1 = [0.49, 0.03, 0.09]; % [s_f, y_f, e_f]
    target_2 = [0.29, 0.09, 0.15]; % [s_f, y_f, e_f]

    tol = 0.01; % convergence tolerance
    max_iter = 100;

    % newton raphson loop
    for iter = 1:max_iter
        [errors, jacobian] = find_errors_and_jacobian(seeds, target_1, ...
            tspan, a, b, c, d, g, x2_crit, d_values);

        % Check convergence
        if all(abs(errors) < tol)
            break;
        end

        dx = -jacobian \ errors;
        % update
        guess = guess + dx;
    end


    % plot everything
    final = [guess(1), guess(2), guess(3), 0.0];
    x_final = rk4(@given_model, final, tspan, a, b, c, d, g, x2_crit);
    
    figure;
    plot(tspan, x_final);
    legend('Sugar', 'Oxygen', 'Yeast', 'Ethanol');
    xlabel('time');
    ylabel('amount');
    title('title');
    grid on;
end


% found nearly entire thing from a MATLAB forum
% https://www.mathworks.com/matlabcentral/answers/460395-runge-kutta-4th-order-method#accepted_answer_373751
function x = rk4(model, x0, tspan, a, b, c, d, g, x2_crit)
    N = length(tspan);
    dt = tspan(2) - tspan(1);
    x = zeros(N, length(x0));
    x(1,:) = x0';
    
    for i = 1:N-1
        k1 = model(tspan(i)     , x(i,:)'            , a, b, c, d, g, x2_crit);
        k2 = model(tspan(i)+dt/2, x(i,:)' + dt/2 * k1, a, b, c, d, g, x2_crit);
        k3 = model(tspan(i)+dt/2, x(i,:)' + dt/2 * k2, a, b, c, d, g, x2_crit);
        k4 = model(tspan(i)+dt  , x(i,:)' + dt * k3  , a, b, c, d, g, x2_crit);
        
        x(i+1,:) = x(i,:) + (dt/6) * (k1 + 2 * k2 + 2 * k3 + k4)';
    end
end


% given formulas from writeup
function solution = given_model(t, x, a, b, c, d, g, x2_crit)
    solution = zeros(4,1);
    solution(1) = -a * x(3);
    solution(2) = -b * x(1) * x(3);
    solution(3) = c * x(1) * x(3) * (x(2) - x2_crit) - d * x(4);
    solution(4) = g * x(1) * x(3);
end


function [error_values, Jacobian_matrix] = find_errors_and_jacobian(seeds, ...
    target_values, tspan, a, b, c, d, g, x2_crit, d_values)

    starting_values = [seeds(1); seeds(2); seeds(3); 0.0];
    
    % stock values
    x = rk4(@given_model, starting_values, tspan, a, b, c, d, g, x2_crit);
    x = x(end, :);
    error_values = [x(1) - target_values(1); x(3) - target_values(2); x(4) - target_values(3)];

    % perturbed values
    Jacobian_matrix = zeros(3,3);
    % calculating the x values with sugar preturbed
    x_s = rk4(@given_model, [seeds(1) + d_values(1); seeds(2); seeds(3); 0.0], ...
        tspan, a, b, c, d, g, x2_crit);
    % assigning to appropriate row in jacobian matrix
    Jacobian_matrix(:, 1) = ([x_s(end, 1); x_s(end, 3); x_s(end, 4)] - [x(1); x(3); x(4)]) / d_values(1);

    % calculating the x values with yeast preturbed
    x_y = rk4(@given_model, [seeds(1); seeds(2) + d_values(2); seeds(3); 0.0], ...
        tspan, a, b, c, d, g, x2_crit);
    % assigning to appropriate row in jacobian matrix
    Jacobian_matrix(:, 2) = ([x_y(end, 1); x_y(end, 3); x_y(end, 4)] - [x(1); x(3); x(4)]) / d_values(2);

    % calculating the x values with ethanol preturbed
    x_e = rk4(@given_model, [seeds(1); seeds(2); seeds(3) + d_values(3); 0.0], ...
        tspan, a, b, c, d, g, x2_crit);
    % assigning to appropriate row in jacobian matrix
    Jacobian_matrix(:, 3) = ([x_e(end, 1); x_e(end, 3); x_e(end, 4)] - [x(1); x(3); x(4)]) / d_values(3);

end


