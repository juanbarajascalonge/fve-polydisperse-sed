%% Example 3 (Section 4.3): N = 2, trapezoidal (roof-shaped) vessel
%  Meshes/roof_<k>.mat, k = 1,...,4  (h = 0.5, 0.25, 0.125, 0.0625).
%  Output: Results/ex3_k<k>_<fve|irp>_t<time>.mat

clear; close all; clc
addpath('src');
levels = 1:4;
scheme = {'fve', 'irp'};

model = default_model('d', [2.9 2.0]*1e-3, 'nrz', 4.6, 'phimax', 0.6, ...
                      'mu0', 100, 'k', [1 0]);
Phi0  = [0.06 0.02];

for k = levels
    mesh = load_mesh(sprintf('Meshes/roof_%d.mat', k));
    for s = 1:2
        tf = 1;
        opt = struct('tf', tf, 'save_times', 1, 'irp', s == 2, ...
                     'name', sprintf('Results/ex3_k%d_%s', k, scheme{s}));
        solve_fve(repmat(Phi0, numel(mesh.g.A), 1), mesh, model, opt);
    end
end
