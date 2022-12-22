function event2txt(event, data_name)

event.t = event.t - event.t(1) + 1;
h       = event.h;
w       = event.w;
n_event = size(event.x, 1);

fileID  = fopen(['data/', data_name, '.txt'], 'w');
fprintf(fileID, '%s,%s,%s,%s,%s,%s\n', 'x', 'y', 't', 'p', 'vx_perp', 'vy_perp');
for i = 1 : n_event
    if event.t(i) < inf
        fprintf(fileID, '%d,%d,%d,%d,%f,%f\n', ...
            event.x(i)-1, event.y(i)-1, int64(event.t(i)), event.p(i), event.vx_perp(i), event.vy_perp(i));
    end
end
fclose(fileID);

end