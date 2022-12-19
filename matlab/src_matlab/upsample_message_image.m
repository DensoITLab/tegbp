function nodes = upsample_message_image(nodes, l)

% nodes{l-1}.m_lr     = (nodes{l-1}.m_lr + repelem(nodes{l}.m_lr, 1, 2, 2)) / 2;
% nodes{l-1}.m_rl     = (nodes{l-1}.m_rl + repelem(nodes{l}.m_rl, 1, 2, 2)) / 2;
% nodes{l-1}.m_ud     = (nodes{l-1}.m_ud + repelem(nodes{l}.m_ud, 1, 2, 2)) / 2;
% nodes{l-1}.m_du     = (nodes{l-1}.m_du + repelem(nodes{l}.m_du, 1, 2, 2)) / 2;

nodes{l-1}.m_lr     = repelem(nodes{l}.m_lr, 1, 2, 2);
nodes{l-1}.m_rl     = repelem(nodes{l}.m_rl, 1, 2, 2);
nodes{l-1}.m_ud     = repelem(nodes{l}.m_ud, 1, 2, 2);
nodes{l-1}.m_du     = repelem(nodes{l}.m_du, 1, 2, 2);

end