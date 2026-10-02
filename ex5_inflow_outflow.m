%% Example 5 (Section 4.5): N = 4, vessel with continuous inflow and outflow
%  Rectangular tank of length 4 and height 1, mouths on the upper 0.3 of the
%  left (inflow) and right (outflow) walls, Poiseuille profiles with mean
%  velocity 2.5, feed Phi_in = (0.03,0.03,0.03,0.03), clear fluid at t = 0,
%  final time t = 4.25 (Fig. 8).
%  Coordinates: x along the tank, y upwards (gravity -y); the paper uses x
%  downwards.  Mesh: 160 x 30 cells, 9.600 triangles, 14.590 control volumes.
%  Output: Results/ex5_t<time>.mat

clear; close all; clc
addpath('src');
mesh_file = 'Meshes/tank_160x30.mat';
L = 4;  H = 1;  nx = 160;  ny = 30;
geom   = struct('x0', 0, 'x1', L, 'yin', [0.7 1.0], 'yout', [0.7 1.0], 'tol', 1e-10*max(L,H));
Ubar   = 2.5;
Phi_in = [0.03 0.03 0.03 0.03];
tf     = 4.25;
model  = default_model('d', [2.9 2.5 2.0 1.3]*1e-3, 'nrz', 4.6, 'phimax', 0.6, ...
                       'mu0', 100, 'k', [0 -1]);

% structured triangulation with alternating diagonals (generated once)
if ~exist(mesh_file, 'file')
    [X, Y] = ndgrid(linspace(0,L,nx+1), linspace(0,H,ny+1));
    Vtri = [X(:).'; Y(:).'];
    [J, I] = ndgrid(1:ny, 1:nx);  I = I(:).';  J = J(:).';
    id = @(i,j) (j-1)*(nx+1) + i;
    a = id(I,J);  b = id(I+1,J);  c = id(I+1,J+1);  d = id(I,J+1);
    alt = mod(I+J, 2) == 0;
    T1 = [a; b; c];  T1(:,~alt) = [a(~alt); b(~alt); d(~alt)];
    T2 = [a; c; d];  T2(:,~alt) = [b(~alt); c(~alt); d(~alt)];
    Ttri = sort(reshape([T1; T2], 3, []), 1);
    save(mesh_file, 'Vtri', 'Ttri');
end
mesh = load_mesh(mesh_file);
bc   = inflow_outflow_bc(mesh, geom, Ubar, Phi_in);

opt = struct('tf', tf, 'save_times', [1 2 3 4 tf], 'bc', bc, 'name', 'Results/ex5');
solve_fve(zeros(numel(mesh.g.A), model.N), mesh, model, opt);
