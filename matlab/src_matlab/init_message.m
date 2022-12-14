function node = init_message(node, idx, idx_4n, act_4n)

%% initialize msg

%% inact"からの"msgをinit
node.m_lr(:, idx(~act_4n(1, :)))    = 0;
node.m_rl(:, idx(~act_4n(2, :)))    = 0;
node.m_du(:, idx(~act_4n(3, :)))    = 0;
node.m_ud(:, idx(~act_4n(4, :)))    = 0;

%% inact　init
node.marg(:, idx_4n(~act_4n))       = 0;
node.m_yx(:, idx_4n(~act_4n))       = 0;

end