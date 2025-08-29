function stent_development()
    L_vec = 5 : 30;
    D_vec = 10 : 50;
    N_total = 1000;

    % Initialize probability matrix
    P = zeros(length(D_vec), length(L_vec));

    % Loop over all (L, D) pairs
    for i = 1:length(L_vec)
        for j = 1:length(D_vec)
            L = L_vec(i);
            D = D_vec(j);
            
            % Generate N_total random S values
            S = 0.02 * D^2 * L + 6 * randn(1, N_total);
            
            % Compute effectiveness E for all generated S values
            E = 24 * log(L * D^2) + 0.18 * D^2 * L - 9.5 * (S + 4);
            
            % Compute probability of successful treatment (E >= 100)
            P(j, i) = sum(E >= 100) / N_total;
        end
    end

    [L_plot, D_plot] = meshgrid(L_vec, D_vec);
    figure;
    surf(L_plot, D_plot, P);
    xlabel('Stent Lifetime L');
    ylabel('Drug Dosage D');
    zlabel('Success Probability P');
    title('Stent Effectiveness');
    colorbar;
    grid on;
end
