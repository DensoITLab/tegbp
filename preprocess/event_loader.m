function event = event_loader(path, padding)

l = load(path);
event.x = int16(l.event.x) + 1;
event.y = int16(l.event.y) + 1;
event.t = l.event.ts;
event.p = l.event.p;
event.h = max(event.y);
event.w = max(event.x);

if padding > 0
    event.x = event.x + padding;
    event.y = event.y + padding;
    event.h = event.h + 2*padding;
    event.w = event.w + 2*padding;
end

end