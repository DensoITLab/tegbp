function msg = smoothness_factor(come, eta, Lam, state, huber)

% messgae passing at factor node
% smoothness term

% 概要：今の状態から近かったらガウス、遠かったらlinearになるようにスケールする(huber)
% 今回のsmoothness factorは差分, hは差分

if huber > 0
    M2      = state' * Lam * state; % linear (Lam prime)
    M       = sqrt(M2);
    if M2 < huber*huber
        scale = 1;
    else
        scale   = 2*huber./M - huber*huber./M2;
    end
else
    scale   = 1;
end
eta     = eta * scale;
Lam     = Lam * scale;

%% calculate message
Lam_aa      = Lam(1:2, 1:2, :);
Lam_ab      = Lam(1:2, 3:4, :);
Lam_bb      = Lam(3:4, 3:4) + reshape(come(3:end), 2, 2);
eta_a       = eta(1:2, :, :);
eta_b       = eta(3:4, :, :) + come(1:2, :);

eta         = eta_a - Lam_ab / Lam_bb * eta_b;
Lam         = Lam_aa - Lam_ab / Lam_bb * Lam_ab';

msg         = [eta; Lam(:)];

end