function [F, t] = droplet_radiation(c, T_ship, T_initial, T_final)
    c = 4.3 * 10^-10;
    T_ship = 290;
    T_initial = 500;
    T_final = 300;

    F = 0 : 0.05 : 1;
    t = zeros(size(F));  % Initialize time array

    % Loop over different values of F
    for i = 1:length(F)
        % Define the function to integrate
        func = @(T) 1 ./ (-c * (1 - F(i)) * T.^4 - c * F(i) * (T.^4 - T_ship^4));

        % Perform numerical integration using the trapezoider function
        t(i) = trapezoider(func, T_initial, T_final, 10);
    end

    % Plot the results
    figure;
    scatter(F, t, 'blue','filled');
    hold on;
    plot(F, t, 'red', 'LineWidth', 2);
    xlabel('View Factor F');
    ylabel('Cooling Time t (seconds)');
    title('Cooling Time vs View Factor');
    grid on;

end

function [A] = trapezoider(f, x_L, x_U, N)
    % following https://byjus.com/maths/trapezoidal-rule/

    % first have to get the step size
    % delta x = (b - a) / N
    delta_x = (x_U - x_L) / N;

    % get all the x values that we need
    x_values = linspace(x_L, x_U, N+1);

    % get all the y values that we need
    y_values = f(x_values);

    % apply the trapezoidal rule
    A = (delta_x/2) * (y_values(1) + 2*sum(y_values(2:end-1)) + y_values(end));
end