function nodes = upsample_message(nodes, y, x, l)

% nodes{l-1}.m_lr     = (nodes{l-1}.m_lr + repelem(nodes{l}.m_lr, 1, 2, 2)) / 2;
% nodes{l-1}.m_rl     = (nodes{l-1}.m_rl + repelem(nodes{l}.m_rl, 1, 2, 2)) / 2;
% nodes{l-1}.m_ud     = (nodes{l-1}.m_ud + repelem(nodes{l}.m_ud, 1, 2, 2)) / 2;
% nodes{l-1}.m_du     = (nodes{l-1}.m_du + repelem(nodes{l}.m_du, 1, 2, 2)) / 2;

x_patch_from    = (-int32(1) : int32(1)) + x;
y_patch_from    = (-int32(1) : int32(1)) + y;

x_patch_to      = (-int32(3) : int32(2)) + 2*x;
y_patch_to      = (-int32(3) : int32(2)) + 2*y;

nodes{l-1}.m_lr(:, y_patch_to, x_patch_to) = ...
    repelem(nodes{l}.m_lr(:, y_patch_from, x_patch_from), 1, 2, 2);
nodes{l-1}.m_rl(:, y_patch_to, x_patch_to) = ...
    repelem(nodes{l}.m_rl(:, y_patch_from, x_patch_from), 1, 2, 2);
nodes{l-1}.m_ud(:, y_patch_to, x_patch_to) = ...
    repelem(nodes{l}.m_ud(:, y_patch_from, x_patch_from), 1, 2, 2);
nodes{l-1}.m_du(:, y_patch_to, x_patch_to) = ...
    repelem(nodes{l}.m_du(:, y_patch_from, x_patch_from), 1, 2, 2);

end