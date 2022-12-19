function mu = belief_vec_to_mu(belief)

eta     = belief(1:2, :);
Lam     = reshape(belief(3:6, :), 2, 2);
mu      = Lam \ eta;

end