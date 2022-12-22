function save_normal_flow_event(event, data_name)

ind_notnan      = ~isnan(event.vx_perp(:, 1));
event.t         = event.t(ind_notnan);
event.x         = event.x(ind_notnan);
event.y         = event.y(ind_notnan);
event.p         = event.p(ind_notnan);
event.vx_perp   = event.vx_perp(ind_notnan);
event.vy_perp   = event.vy_perp(ind_notnan);

event2txt(event, data_name);

end