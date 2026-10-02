%% Example 4 (Section 4.4): N = 4, reflux-classifier geometry
%  Three cases up to t = 4 (Fig. 7):
%    (a) mu0 = 100,  n_l = 4.6
%    (b) mu0 = 100,  (n_1,...,n_4) = (3.1, 3.5, 4.0, 4.6)
%    (c) mu0 = 1000, n_l = 4.6
%  Output: Results/ex4<a|b|c>_t<time>.mat

clear; close all; clc
addpath('src');
mesh  = load_mesh('Meshes/refluxclas_3.mat');
cases = struct('tag', {'a', 'b', 'c'}, 'mu0', {100, 100, 1000}, ...
               'nrz', {4.6, [3.1 3.5 4.0 4.6], 4.6});
Phi0  = [0.05 0.05 0.05 0.05];
tf    = 4;

for c = cases
    model = default_model('d', [2.9 2.5 2.0 1.3]*1e-3, 'nrz', c.nrz, ...
                          'phimax', 0.6, 'mu0', c.mu0, 'k', [0 -1]);
    opt = struct('tf', tf, 'save_times', 1:tf, 'name', ['Results/ex4' c.tag]);
    solve_fve(repmat(Phi0, numel(mesh.g.A), 1), mesh, model, opt);
end
