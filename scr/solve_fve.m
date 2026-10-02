function out = solve_fve(Phi0, mesh, model, opt)
%SOLVE_FVE  Coupled FVE scheme of Algorithm 3.1.
%
%   out = solve_fve(Phi0, mesh, model, opt)
%
%   Phi0  : NC x N initial cell averages     mesh : LOAD_MESH / BUILD_MESH
%   model : DEFAULT_MODEL
%   opt fields (defaults in brackets):
%     tf           final time
%     cfl          CFL number of Algorithm 3.1                    
%     rk           2 = SSPRK2 (paper), 1 = forward Euler           
%     irp          true = IRP-FVE (scaling limiters), false = FVE  
%     save_times   times at which a snapshot is written           
%     name         snapshots go to <name>_t<time>.mat ('' = none)  
%     bc           struct from INFLOW_OUTFLOW_BC (Example 5)       
%     source       @(x,y,t) NC x N source of the transport equations  (Example 1)
%     force        @(x,y,t) NE x 2 extra force of the Stokes problem (Example 1)
%     verbose      print every 'verbose' steps                     


cfl  = getopt(opt, 'cfl', 0.45);
rk   = getopt(opt, 'rk', 2);
irp  = getopt(opt, 'irp', true);
name = getopt(opt, 'name', '');
bc   = getopt(opt, 'bc', []);
src  = getopt(opt, 'source', []);
frc  = getopt(opt, 'force', []);
every = getopt(opt, 'verbose', 100);
tf   = opt.tf;
ts   = sort(getopt(opt, 'save_times', tf));
ts   = ts(ts > 0 & ts <= tf + 1e-12);
if isempty(ts) || abs(ts(end) - tf) > 1e-12, ts = [ts(:).', tf]; end

if isempty(bc)
    faceTag = [];  Phi_in = [];
else
    faceTag = bc.faceTag;  Phi_in = bc.Phi_in;
end
S  = stokes_setup(mesh, bc);
FV = fv_setup(mesh, model, faceTag);
x  = FV.xC;  y = FV.yC;                     % dual cell centres = primal edge midpoints

rhs = @(Phi, t) stage(Phi, t, S, FV, model, irp, Phi_in, src, frc);
g = mesh.g;  gt = mesh.gt;  gt.elem2dof = mesh.elem2dof;   % stored with every snapshot

u = Phi0;  t = 0;  n = 0;  is = 1;
minphi = min(u, [], 1);  maxphi = max(sum(u, 2));
fprintf('%10s %11s %13s %11s %11s\n', 't', 'dt', 'min phi_l', 'max phi', 'max|div|');
tic0 = tic;
while t < tf - 1e-12
    [R1, alpha, divmax] = rhs(u, t);
    dt = cfl*FV.dt_geo/max(alpha, eps);
    dt = min(dt, ts(is) - t);
    if rk == 1
        u = u + dt*R1;
    else
        u1 = u + dt*R1;
        u  = 0.5*(u + u1 + dt*rhs(u1, t + dt));
    end
    t = t + dt;  n = n + 1;
    if abs(t - ts(is)) < 1e-12, t = ts(is); end
    minphi = min(minphi, min(u, [], 1));
    maxphi = max(maxphi, max(sum(u, 2)));
    if mod(n, every) == 0 || abs(t - ts(is)) < 1e-12
        fprintf('%10.4f %11.4e %13.4e %11.6f %11.2e\n', t, dt, min(u(:)), max(sum(u,2)), divmax);
    end
    if abs(t - ts(is)) < 1e-12
        if ~isempty(name)
            [qx, qy, p] = stokes_solve(S, model, u, evalf(frc, x, y, t));
            file = sprintf('%s_t%s.mat', name, strrep(num2str(ts(is)), '.', 'p'));
            save(file, 'u', 'qx', 'qy', 'p', 't', 'minphi', 'maxphi', 'g', 'gt');
            fprintf('  -> %s\n', file);
        end
        is = min(is + 1, numel(ts));
    end
end
fprintf('%d steps, %.1f s\n', n, toc(tic0));
[qx, qy, p] = stokes_solve(S, model, u, evalf(frc, x, y, t));
out = struct('u', u, 'qx', qx, 'qy', qy, 'p', p, 't', t, ...
             'minphi', minphi, 'maxphi', maxphi, 'nsteps', n);
end

function [R, alpha, divmax] = stage(Phi, t, S, FV, model, irp, Phi_in, src, frc)
[ux, uy, ~, info] = stokes_solve(S, model, Phi, evalf(frc, FV.xC, FV.yC, t));
[R, alpha] = fv_residual(Phi, ux, uy, FV, model, irp, Phi_in);
if ~isempty(src), R = R + src(FV.xC, FV.yC, t); end
divmax = info.divmax;
end

function F = evalf(f, x, y, t)
if isempty(f), F = []; else, F = f(x, y, t); end
end

function v = getopt(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
