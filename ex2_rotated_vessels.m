%% Example 2 (Section 4.2): N = 2, effect of the limiters
%  Rectangular vessels rotated by theta = 0, 30, 45, 60 degrees; meshes
%  Meshes/inclined_<theta>_<k>.mat, k = 1,...,4  (h = 0.5, 0.25, 0.125, 0.0625).
%  Output: Results/ex2_th<theta>_k<k>_<fve|irp>_t<time>.mat

clear; close all; clc
addpath('src');
angles = [0 30 45 60];
levels = 1:4;
scheme = {'fve', 'irp'};                       % without / with scaling limiters

model = default_model('d', [2.9 2.0]*1e-3, 'nrz', 4.6, 'phimax', 0.6, ...
                      'mu0', 1000, 'k', [1 0]);
Phi0  = [0.06 0.02];

for th = angles
    for k = levels
        mesh = load_mesh(sprintf('Meshes/inclined_%d_%d.mat', th, k));
        for s = 1:2
            tf = 1.5;
            opt = struct('tf', tf, 'save_times', [0.5 1.5 3], 'irp', s == 2, ...
                         'name', sprintf('Results/ex2_th%d_k%d_%s', th, k, scheme{s}));
            solve_fve(repmat(Phi0, numel(mesh.g.A), 1), mesh, model, opt);
        end
    end
end
