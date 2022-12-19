%% jnagata 2022/12/08
% cpp実装のためにMATLABでminimumな実装する。

%%
clear
close all

addpath src_matlab outer/flow-code-matlab/

%% data load
load("data/bricks/event.mat");
event.t = event.t - event.t(1) + 1;
h       = event.h;
w       = event.w;
n_event = size(event.x, 1);

%% parameter
n_scale         = 1;
sig2_p          = 0.01;
sig2_r          = 9;
sig2_t          = 100;
n_batch         = 100;
huber           = 0;
dt_act          = 50000;
n_iter          = 1;

flg_imgvis      = true;

%% node
nodes           = define_nodes(h, w);

%% prepare
Sig_m           = [sig2_r, 0; 0, sig2_t];
Lam_m           = inv(Sig_m);
Sig_p           = eye(2) * sig2_p;
J               = [-1, 0, 1, 0; 0, -1, 0, 1];
Lam_p           = J' / Sig_p * J;
pad_size        = 16;

%% keep
flow_gt_        = zeros(2, h, w);
flow_est_       = zeros(2, h, w);
flow_norm_      = zeros(2, h, w);
t_act_          = ones(h, w) * -inf;

indices         = zeros(2, n_batch, 'int32');
v_norms         = nan(2, n_batch);
timestamps      = nan(1, n_batch);
cnt             = 0;

%% frame
dt_frame        = 50000;
t_frame_max     = floor(max(event.t) / dt_frame) * dt_frame;
t_frame_list    = dt_frame : dt_frame : t_frame_max;
i_frame         = 1;
n_frame         = length(t_frame_list);

%% display
if flg_imgvis
    fig2    = figure;
    set(fig2, 'Position', [1526 14 994 807])
    ax2     = cell(1, 6);
    im2     = cell(1, 6);
    titles  = {'norm', 'Ours'};
    for i = 1 : 2
        ax2{i} = subaxis(1, 2, i, 'S', 0, 'M', 0, 'P', 0.002);
        im2{i} = imagesc(zeros(h-2*pad_size, w-2*pad_size, 'uint8'));
        title(titles{i})
        axis off
        axis image
        clim([0 255])

        set(ax2{i}, 'FontSize', 20)
    end
end

%% event loop
disp('start event loop')
start_frame = tic;
for e = 1 : n_event
    %     if flg_profile
    %         profile viewer
    %         profile clear
    %         profile on
    %     end
    xi          = event.x(e);
    yi          = event.y(e);
    ti          = event.t(e);
    qi          = event.p(e);
    vx_perpi    = event.vx_perp(e);
    vy_perpi    = event.vy_perp(e);
    
    %% store normal flow for comparison
    flow_norm_(1, yi, xi) = vx_perpi;
    flow_norm_(2, yi, xi) = vy_perpi;

    %% visualization timing
    if ti > t_frame_list(i_frame) || e == n_event
        time_frame = toc(start_frame);
        start_frame = tic;
        if rem(i_frame, 1) == 0
            fprintf('frame [%d / %d] %.2f s\n', i_frame, n_frame, time_frame);
        end
        t_frame         =  t_frame_list(i_frame);

        %% permute & unpadding for visualization
        flow_est    = permute(nodes.marg(1:2, pad_size+1:end-pad_size, pad_size+1:end-pad_size), [2, 3, 1]);
        flow_norm   = permute(flow_norm_(:, pad_size+1:end-pad_size, pad_size+1:end-pad_size), [2, 3, 1]);
        t_act       = nodes.t_act(pad_size+1:end-pad_size, pad_size+1:end-pad_size);
        mask        = t_act > t_frame - dt_act;
        flow_est(~cat(3, mask, mask))   = 0;
        flow_norm(~cat(3, mask, mask))  = 0;
        I_est       = flowToColor(flow_est, sqrt(2));
        I_norm      = flowToColor(flow_norm, sqrt(2));
        
        %% visualization
        im2{1}.CData   = I_est;
        im2{2}.CData   = I_norm;
        drawnow

        %% next step
        i_frame     = i_frame + 1;
        if e == n_event || i_frame > n_frame
            break
        end
    end

    %% store in buffer
    nodes.t_act(yi, xi)     = ti;
    cnt                     = cnt + 1;
    v_norms(:, cnt)         = [vx_perpi; vy_perpi];
    indices(:, cnt)         = [yi; xi];
    timestamps(:, cnt)      = ti;

    %% batch分溜まったら
    if cnt >= n_batch 
        %% ここを並列化する？
        for i = 1 : n_batch
            v_perp      = v_norms(:, i);
            y           = indices(1, i);
            x           = indices(2, i);
            t           = timestamps(1, i);
     
            %% data factor
            obs_msg     = calc_data_term(v_perp, Lam_m, huber);
            
            %% setting data factor
            nodes.m_yx(:, y, x) = obs_msg;

            %% message passing
            nodes       = message_passing_event(nodes, y, x, t, dt_act, Lam_p, huber);
            [nodes, ~]  = update_state(nodes, y, x, t, dt_act);
        end

        %% new batch
        cnt = 0;
    end
end
