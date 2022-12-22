function [v_norm, inlier_ratio, n_inlier] = planefit_flow(base, target, n_trial, th_t_inlier)

n_point = size(base, 1);

inlier  = ones(n_point, 1, 'logical');
for n = 1 : n_trial
    pln     = pinv(base(inlier, :)) * target(inlier);
    a       = pln(1);
    b       = pln(2);
    z       = sqrt(a*a + b*b);
    
    d_target        = abs(base * pln - target);
    if th_t_inlier > 0
        inlier          = d_target < th_t_inlier;
    else
        inlier          = d_target < z/2;
    end
end
inlier_ratio    = nnz(inlier) / n_point;
n_inlier        = nnz(inlier);
v_norm          = [a; b] ./ (z*z);

end