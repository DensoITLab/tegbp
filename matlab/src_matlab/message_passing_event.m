function nodes = message_passing_event(nodes, y, x, t, dt_act, Lam_p, huber)

dirc            = [0, 1; 0, -1; 1, 0; -1, 0];

[nodes, belief] = update_state(nodes, y, x, t, dt_act); % collect and update state
nodes           = send_message_4connect(nodes, y, x, belief, t, dt_act, Lam_p, huber); % send messaages

for d = 1 : 4
    y_to = y + dirc(d, 1);
    x_to = x + dirc(d, 2);
    if t - nodes.t_act(y_to, x_to) < dt_act % only active neighbor nodes
        [nodes, belief] = update_state(nodes, y_to, x_to, t, dt_act);
        nodes           = send_message_4connect(nodes, y_to, x_to, belief, t, dt_act, Lam_p, huber);
    end 
end

end