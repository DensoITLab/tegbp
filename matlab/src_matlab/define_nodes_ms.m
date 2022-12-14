function nodes = define_nodes_ms(h, w, n_scale)

nodes           = cell(1, n_scale);

for l = 1 : n_scale
    hh      = h / 2^(l-1);
    ww      = w / 2^(l-1);
    
    msg     = zeros(6, hh, ww);
    
    nodes{l}.m_lr  = msg; % 左からのメッセージ
    nodes{l}.m_rl  = msg; % 右からのメッセージ
    nodes{l}.m_ud  = msg; % 上からのメッセージ
    nodes{l}.m_du  = msg; % 下からのメッセージ
    nodes{l}.m_yx  = msg; % 観測からのメッセージ
    nodes{l}.marg  = msg; % 周辺事後分布
    nodes{l}.t_act = ones(hh, ww)*-inf;
end

end