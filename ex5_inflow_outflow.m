%% Example 5 (Section 4.5): N = 4, vessel with continuous inflow and outflow
%  Output: Results/ex5_t<time>.mat

clear; close all; clc
addpath('src');
mesh_file = 'Meshes/tank_160x30.mat';
mesh = load_mesh(mesh_file);

L = 4;  H = 1;  nx = 160;  ny = 30;
geom   = struct('x0', 0, 'x1', L, 'yin', [0.7 1.0], 'yout', [0.7 1.0], 'tol', 1e-10*max(L,H));
Ubar   = 2.5;
Phi_in = [0.03 0.03 0.03 0.03];
tf     = 4.25;
model  = default_model('d', [2.9 2.5 2.0 1.3]*1e-3, 'nrz', 4.6, 'phimax', 0.6, ...
                       'mu0', 100, 'k', [0 -1]);

bc   = inflow_outflow_bc(mesh, geom, Ubar, Phi_in);

opt = struct('tf', tf, 'save_times', tf, 'bc', bc, 'name', 'Results/ex5');
solve_fve(zeros(numel(mesh.g.A), model.N), mesh, model, opt);
