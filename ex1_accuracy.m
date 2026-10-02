%% Example 1 (Section 4.1): N = 3, numerical order of accuracy
%  Manufactured solution (see src/manufactured_ex1.m) on the unit square,
%  quasi-uniform meshes Meshes/square<n>.mat, h = 1/n = 2^-l, l = 2,...,7,
%  final time t = 0.01.  Output: Results/ex1_n<n>_t0p01.mat and the errors
%  and convergence rates of Table 1.

clear; close all; clc
addpath('src');
levels = 2:7;
tf     = 0.01;

model = default_model('delta', [1 0.8 0.6], 'nrz', 4.6, 'phimax', 1, ...
                      'mu0', 1, 'k', [1 0]);
ex = manufactured_ex1(model);

%% Simulation
for l = levels
    n    = 2^l;
    mesh = load_mesh(sprintf('Meshes/square%d.mat', n));
    Phi0 = cell_average(@(x,y) ex.phi(x,y,0), mesh);
    opt  = struct('tf', tf, 'name', sprintf('Results/ex1_n%d', n), ...
                  'source', ex.source, 'force', ex.force);
    solve_fve(Phi0, mesh, model, opt);
end

%% Errors and convergence rates (Table 1)
%  e(Phi)  = sum_l ||phi_{l,h} - phi_{l,ex}||_{L1}  (exact cell averages)
%  e0(u)   = ||u_h - u_ex||_{0},  e1(u) = ||u_h - u_ex||_{1,h},  e0(p) = ||p_h - p_ex||_{0}
%  (velocity and pressure: one-point rule at the barycentres of the triangles)
nl = numel(levels);
h  = zeros(nl,1);  E = zeros(nl,4);
for i = 1:nl
    n = 2^levels(i);  h(i) = 1/n;
    R = load(sprintf('Results/ex1_n%d_t%s.mat', n, strrep(num2str(tf),'.','p')));
    g = R.g;  gt = R.gt;  e2d = gt.elem2dof;
    xb = gt.B(1,:).';  yb = gt.B(2,:).';  aT = gt.A(:);

    % concentrations
    uex    = cell_average(@(x,y) ex.phi(x,y,tf), struct('g', g));
    E(i,1) = sum(g.A(:).*sum(abs(R.u - uex), 2));

    % velocity in L2: CR function at the barycentre = mean of the 3 dofs
    uh     = [mean(R.qx(e2d),2), mean(R.qy(e2d),2)];
    E(i,2) = sqrt(sum(aT.*sum((uh - ex.u(xb,yb)).^2, 2)));

    % velocity in the broken H1 seminorm: grad psi_i = -2 grad lambda_i
    X1 = gt.V(:,gt.T(1,:)).';  X2 = gt.V(:,gt.T(2,:)).';  X3 = gt.V(:,gt.T(3,:)).';
    lx = [X3(:,1)-X2(:,1), X1(:,1)-X3(:,1), X2(:,1)-X1(:,1)];   % edge opposite vertex i
    ly = [X3(:,2)-X2(:,2), X1(:,2)-X3(:,2), X2(:,2)-X1(:,2)];
    sa = 0.5*(lx(:,2).*ly(:,3) - lx(:,3).*ly(:,2));              % signed area
    Dx = ly./sa;  Dy = -lx./sa;
    gh = [sum(R.qx(e2d).*Dx,2), sum(R.qx(e2d).*Dy,2), ...
          sum(R.qy(e2d).*Dx,2), sum(R.qy(e2d).*Dy,2)];           % u1_x u1_y u2_x u2_y
    E(i,3) = sqrt(sum(aT.*sum((gh - ex.gradu(xb,yb)).^2, 2)));

    % pressure
    E(i,4) = sqrt(sum(aT.*(R.p - ex.p(xb,yb)).^2));
end

r = zeros(nl,4);                                  % h_{i-1}/h_i = 2
r(2:end,:) = log2(E(1:end-1,:)./E(2:end,:));

fprintf('\n%9s %10s %7s %10s %7s %10s %7s %10s %7s\n', 'h', 'e(Phi)', 'r(Phi)', ...
        'e0(u)', 'r0(u)', 'e1(u)', 'r1(u)', 'e0(p)', 'r0(p)');
for i = 1:nl
    fprintf('%9.2e %10.2e %7.3f %10.2e %7.3f %10.2e %7.3f %10.2e %7.3f\n', h(i), ...
            E(i,1), r(i,1), E(i,2), r(i,2), E(i,3), r(i,3), E(i,4), r(i,4));
end