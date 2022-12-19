function [nodes, belief] = update_state(nodes, y, x, t, dt_act)

belief = nodes.m_yx(:, y, x);

if t - nodes.t_act(y+1, x) < dt_act % active
    belief = belief + nodes.m_du(:, y, x);
end
if t - nodes.t_act(y-1, x) < dt_act % active
    belief = belief + nodes.m_ud(:, y, x);
end
if t - nodes.t_act(y, x-1) < dt_act % active
    belief = belief + nodes.m_lr(:, y, x);
end
if t - nodes.t_act(y, x+1) < dt_act % active
    belief = belief + nodes.m_rl(:, y, x);
end

mu = belief_vec_to_mu(belief);

if isnan(mu(1))
    mu = zeros(2, 1);
end
nodes.marg(1:2, y, x) = mu;
nodes.marg(3:6, y, x) = belief(3:6, :);

end