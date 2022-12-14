clear
close all

%%
dataset = 'toy';
scene   = 'ver0';

dt      = 5000;
dt_act  = 50000;

h       = 128;
w       = 128;

[cols, rows]    = meshgrid(1:w, 1:h);
t               = (cols + rows) * 1000;
p               = ones(h, w);
vx              = (cols >= rows) * 1;
vy              = (cols < rows) * 1;

observation     = [cols(:), rows(:), t(:), p(:), vx(:), vy(:)];

T = table(cols(:)-1, rows(:)-1, t(:), p(:), vx(:), vy(:), ...
    'VariableNames', {'x', 'y', 't', 'p', 'vx', 'vy'});
T = sortrows(T, 't');

writetable(T, fullfile('data', dataset, [scene, '.csv']), 'WriteVariableNames', false)