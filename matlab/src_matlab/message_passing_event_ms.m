function nodes = message_passing_event_ms(nodes, y, x, t, dt_act, Lam_p, huber, n_scale)

dirc            = [0, 1; 0, -1; 1, 0; -1, 0];

for l_ = 1 : n_scale
    l                   = n_scale + 1 - l_;
    
    xx                  = idivide(x, 2^(l-1), 'ceil');
    yy                  = idivide(y, 2^(l-1), 'ceil');

    [nodes{l}, belief]  = update_state(nodes{l}, yy, xx, t, dt_act); % collect and update state
    nodes{l}            = send_message_4connect(nodes{l}, yy, xx, belief, t, dt_act, Lam_p, huber); % send messaages

    for d = 1 : 4
        y_to    = yy + dirc(d, 1);
        x_to    = xx + dirc(d, 2);
        if t - nodes{l}.t_act(y_to, x_to) < dt_act             % only active neighbor nodes
            [nodes{l}, belief]  = update_state(nodes{l}, y_to, x_to, t, dt_act);
            nodes{l}            = send_message_4connect(nodes{l}, y_to, x_to, belief, t, dt_act, Lam_p, huber);
        end
    end
    % upsample
    if l > 1
        nodes = upsample_message(nodes, yy, xx, l);
    end
end

[nodes{l}, ~]  = update_state(nodes{l}, y, x, t, dt_act);

end