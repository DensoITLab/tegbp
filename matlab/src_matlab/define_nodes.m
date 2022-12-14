function nodes = define_nodes(h, w)

msg  = zeros(6, h, w);

nodes.m_lr  = msg; % 左からのメッセージ
nodes.m_rl  = msg; % 右からのメッセージ
nodes.m_ud  = msg; % 上からのメッセージ
nodes.m_du  = msg; % 下からのメッセージ
nodes.m_yx  = msg; % 観測からのメッセージ
nodes.marg  = msg; % 周辺事後分布
nodes.t_act = ones(h, w)*-inf;

end