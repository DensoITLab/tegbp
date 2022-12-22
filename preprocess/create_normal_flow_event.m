% 2022-12-22 jnagata
% calculate normal flow by planar fitting and dump to txt

close all
clearvars

addpath matlab/outer/flow-code-matlab/

%% config (ymlにする?)
win_x       = 2;        % spartial window size
win_t       = 40000;    % temporal window size
t_rf        = 40000;    % refractory filter
n_trial     = 2;        % number of trial for flane fitting
n_th_inlier = 5;        % event number threshold foutlier rejection in plane fitting
t_th_inlier = -1;       % temporal threshold foutlier rejection in plane fitting
th_inlier   = 0.5;      % threshold inlier ratio
dt_scale    = 5000;     % dt scale

data_dir    = '/home/jnagata/petaco/data/DNW/221108_EVK2_150mm_f8/500lux';
data_name   = 'recording_2022-11-08_14-02-40';
save_name   = 'qr000';

%% load data
event_path  = fullfile(data_dir, [data_name, '_cd.mat']);
event       = event_loader(event_path, 16);
h           = event.h;
w           = event.w;
n_event     = length(event.x);

%% prepare
saen        = ones(h, w) * -inf;
saep        = ones(h, w) * -inf;
sae         = ones(h, w) * -inf;
t_fit       = ones(h, w) * -inf;
patchq      = int16(-win_x : win_x);
[bx, by]    = meshgrid(-win_x : win_x);

%% visualize
i_frame     = 1;
batch       = 10^5;
flow_norm   = zeros(h, w, 2);
max_dist    = 100;
max_flow    = sqrt(2);
fig         = figure('Position', [3145 737 1616 420]);
ax          = cell(1, 2);
im          = cell(1, 2);
tl          = tiledlayout(1, 2, 'TileSpacing', 'compact', 'padding', 'compact');
for i = 1 : 2
    ax{i}       = nexttile;
    im{i}       = imagesc(zeros(h, w, 'uint8'));
    axis image
    axis off
end

%% save normal flow
event.vx_perp = nan(n_event, 1);
event.vy_perp = nan(n_event, 1);

%% event loop
disp('start event loop')
t_start     = tic;
ts_start    = event.t(1);
for e = 1 : n_event
    if rem(e, 100000) == 0 % 1M events
        t_process   = toc(t_start);
        t_actual    = event.t(e);
        fprintf('%dK events: processing time %.2f s, actual time %.2f s\n', batch/1000, t_process, t_actual*10^-6);
        
        mask = t_fit > (ti - win_t);

        I_sae = uint8((sae - (ti - win_t)) / win_t * 255); % 正規化
        I_sae(I_sae < 0)    = 0;
        im{1}.CData         = I_sae;

        flow_norm(~cat(3, mask, mask)) = 0; 
        I_norm          = flowToColor(flow_norm, max_flow);
        im{2}.CData     = I_norm;

        drawnow

        %% next slide
        t_start     = tic;
    end

    xi = event.x(e);
    yi = event.y(e);
    ti = event.t(e);
    qi = event.p(e);
    
    if (ti - sae(yi, xi)) < t_rf
        continue
    end

    if qi > 0
        saep(yi, xi) = ti;
    else
        saen(yi, xi) = ti;
    end
    sae(yi, xi) = ti;
    
    if ti < win_t % pass until time window
        continue
    end
    Xq      = xi + patchq;
    Yq      = yi + patchq;
    if qi > 0
        sae_ref = (saep(Yq, Xq) - ti);
    else
        sae_ref = (saen(Yq, Xq) - ti);
    end

    pf_valid    = (sae_ref > - win_t);
    n_point     = nnz(pf_valid);
    
    if n_point >= n_th_inlier
        B       = [bx(pf_valid), by(pf_valid), ones(n_point, 1)];
        f       = sae_ref(pf_valid);
       
        [v_norm, inlier_ratio, n_inlier] = planefit_flow(B, f, n_trial, t_th_inlier);
        v_norm  = v_norm * dt_scale;
        
        if isnan(v_norm(1))
            continue
        end
        if inlier_ratio < th_inlier || n_inlier < n_th_inlier
            continue
        end
        % 大きすぎるのは除く
        if (sqrt(sum(v_norm.^2)) / dt_scale * 10^6 / 30) > max_dist
            continue
        end 

        flow_norm(yi, xi, 1)    = v_norm(1);
        flow_norm(yi, xi, 2)    = v_norm(2);
        event.vx_perp(e)        = v_norm(1);
        event.vy_perp(e)        = v_norm(2);
        t_fit(yi, xi)           = ti;
    end

end

save_normal_flow_event(event, save_name)


