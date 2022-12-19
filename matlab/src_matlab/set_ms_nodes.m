function nodes = set_ms_nodes(nodes, y, x, t, obs_msg, dt_act)

dirc        = [0, 0; 0, -1; -1, 0; -1, -1];

n_scale     = size(nodes, 2);

x_to        = x;
y_to        = y;

%% data項
nodes{1}.m_yx(:, y, x) = obs_msg;

%% 上の解像度
for l = 2 : n_scale
    % 2x2 patch内
    x_from              = x_to;
    y_from              = y_to;
    x_to                = idivide(x_from, 2, 'ceil');
    y_to                = idivide(y_from, 2, 'ceil');
    
    nodes{l}.t_act(y_to, x_to) = nodes{l-1}.t_act(y_from, x_from);
    
    obs_node_act        = zeros(6, 1);
    for k = 1 : 4
        x_from_patch    = x_to * 2 + dirc(k, 1);
        y_from_patch    = y_to * 2 + dirc(k, 2);
        
        if nodes{l-1}.t_act(y_from_patch, x_from_patch) > (t-dt_act)
            obs_node_act = obs_node_act + nodes{l-1}.m_yx(:, y_from_patch, x_from_patch);
        end
    end
    nodes{l}.m_yx(:, y_to, x_to) = obs_node_act;
end

end