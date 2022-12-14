function obs_msg = calc_data_term(v_perp, Lam, huber)
%% 
% Input:
% v_perp    [2 x 1] : normal flow 
% Lam_m     [2 x 2] : measurment 
% Output:
% obs_msg   [6 x 1] : 

theta   = atand(v_perp(2) ./ v_perp(1));
rot_mat = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
Lam_teg = rot_mat * Lam * rot_mat';
eta     = Lam_teg * v_perp;
obs_msg = [eta; Lam_teg(:)];

%% huber
if huber > 0
    M2      = v_perp' * Lam * v_perp; % linear (Lam prime)
    M       = sqrt(M2);
    if M2 < huber*huber
        scale   = 1;
    else
        scale   = 2*huber./M - huber*huber./M2;
    end
else
    scale   = 1;
end
obs_msg     = obs_msg .* scale;

end