function nodes = send_message_4connect(nodes, y, x, belief, t, dt_act, Lam_p, huber)

eta_p       = zeros(4, 1); % z = 0だから    
state_from  = nodes.marg(1:2, y, x);

%% L → R
if t - nodes.t_act(y, x+1) < dt_act % active
    msg_v = belief - nodes.m_rl(:, y, x);
else
    msg_v = belief;
end
state_to    = nodes.marg(1:2, y, x+1);
state       = [state_to; state_from];
msg_p       = smoothness_factor(msg_v, eta_p, Lam_p, state, huber);
nodes.m_lr(:, y, x+1) = msg_p;

%% R → L
if t - nodes.t_act(y, x-1) < dt_act % active
    msg_v = belief - nodes.m_lr(:, y, x);
else
    msg_v = belief;
end
state_to    = nodes.marg(1:2, y, x-1);
state       = [state_to; state_from];
msg_p       = smoothness_factor(msg_v, eta_p, Lam_p, state, huber);
nodes.m_rl(:, y, x-1) = msg_p;

%% U → D
if t - nodes.t_act(y+1, x) < dt_act % active
    msg_v = belief - nodes.m_du(:, y, x);
else
    msg_v = belief;
end
state_to    = nodes.marg(1:2, y+1, x);
state       = [state_to; state_from];
msg_p       = smoothness_factor(msg_v, eta_p, Lam_p, state, huber);
nodes.m_ud(:, y+1, x) = msg_p;

%% D → U
if t - nodes.t_act(y-1, x) < dt_act % active
    msg_v = belief - nodes.m_ud(:, y, x);
else
    msg_v = belief;
end
state_to    = nodes.marg(1:2, y-1, x);
state       = [state_to; state_from];
msg_p       = smoothness_factor(msg_v, eta_p, Lam_p, state, huber);
nodes.m_du(:, y-1, x) = msg_p;


end