function plot_result(res, field, varargin)
%PLOT_RESULT  Plot one field of a snapshot written by the example scripts.
%
%   plot_result(res, field)
%   plot_result(res, field, 'Name', value, ...)
%
%   res    file name of a snapshot (e.g. 'Results/ex3_k4_irp_t3.mat') or the
%          struct obtained with load(...)
%   field  'phi1', ..., 'phiN'  concentration of species l   (dual mesh)
%          'phi'                total concentration            (dual mesh)
%          'ux', 'uy', 'umag'   velocity components / |u|      (primal mesh,
%                               value at the barycentre of each triangle)
%          'p'                  pressure                       (primal mesh)
%   Options:
%     'view'     [az el] passed to view()            [0 90]; use [90 90] when
%                                                     x points downwards
%     'rotate'   rotation of the domain in degrees    0
%     'clim'     colour limits                        automatic
%     'colorbar' colorbar location, or 'none'         'eastoutside'
%     'quiver'   draw every k-th velocity arrow, 0 = none   0
%     'export'   file name (.eps, .png, ...) to save the figure   ''
%
%   Examples
%     plot_result('Results/ex3_k4_irp_t3.mat', 'phi', 'view', [90 90])
%     plot_result('Results/ex4a_t4.mat', 'umag', 'rotate', 90, 'view', [90 90])
%     plot_result('Results/ex5_t4p25.mat', 'phi1', 'colorbar', 'southoutside', ...
%                 'export', 'Figures/ex5_phi1.eps')

opt = struct('view', [0 90], 'rotate', 0, 'clim', [], 'colorbar', 'eastoutside', ...
             'quiver', 0, 'export', '');
for i = 1:2:numel(varargin), opt.(lower(varargin{i})) = varargin{i+1}; end
if ischar(res), res = load(res); end
g = res.g;  gt = res.gt;

% cellwise values and the mesh they live on
if strncmp(field, 'phi', 3)
    if strcmp(field, 'phi'), U = sum(res.u, 2); else, U = res.u(:, str2double(field(4:end))); end
    V = double(g.V);  T = double(g.T);
else
    V = double(gt.V);  T = double(gt.T);
    e2d = elem2dof(gt);
    ux = mean(res.qx(e2d), 2);  uy = mean(res.qy(e2d), 2);
    switch field
        case 'ux',   U = ux;
        case 'uy',   U = uy;
        case 'umag', U = sqrt(ux.^2 + uy.^2);
        case 'p',    U = res.p;
        otherwise,   error('plot_result:field', 'Unknown field ''%s''.', field);
    end
end
R = [cosd(opt.rotate) -sind(opt.rotate); sind(opt.rotate) cosd(opt.rotate)];
V = R*V;

fig = figure('Color', 'w');
F = T.';  F(F <= 0) = NaN;                       % dual cells: T(4,:) = -1 on triangles
patch('Faces', F, 'Vertices', V.', 'FaceVertexCData', U(:), ...
      'FaceColor', 'flat', 'EdgeColor', 'none');
hold on

% domain boundary: polygon sides that belong to a single cell
nv = sum(T > 0, 1);  E = zeros(0,2);
for k = 1:size(T,1)
    c = find(nv >= k);
    E = [E; T(k,c).', T(sub2ind(size(T), mod(k, nv(c)) + 1, c)).'];   %#ok<AGROW>
end
[ue, ~, ic] = unique(sort(E,2), 'rows');
be = ue(accumarray(ic, 1) == 1, :);
X = [V(1,be(:,1)); V(1,be(:,2)); nan(1,size(be,1))];
Y = [V(2,be(:,1)); V(2,be(:,2)); nan(1,size(be,1))];
plot(X(:), Y(:), 'k-', 'LineWidth', 1);

if opt.quiver > 0
    if strncmp(field, 'phi', 3)
        e2d = elem2dof(gt);  ux = mean(res.qx(e2d), 2);  uy = mean(res.qy(e2d), 2);
    end
    B = R*double(gt.B);  W = R*[ux.'; uy.'];  W = W/max(sqrt(sum(W.^2, 1)));
    k = 1:opt.quiver:size(B,2);
    quiver(B(1,k), B(2,k), W(1,k), W(2,k), 1, 'w');
end
hold off

axis image off
view(opt.view);
colormap(turbo);
if ~isempty(opt.clim), caxis(opt.clim); end
if ~strcmp(opt.colorbar, 'none'), colorbar(opt.colorbar); end

if ~isempty(opt.export)
    if exist('exportgraphics', 'file') || exist('exportgraphics', 'builtin')
        exportgraphics(fig, opt.export, 'Resolution', 600, 'ContentType', 'image');
    else
        print(fig, opt.export, '-depsc', '-r600');
    end
end
end

%% ------------------------------------------------------------------------
function e2d = elem2dof(gt)
% CR dofs of each triangle (local edge i opposite vertex i); recomputed if the
% snapshot does not store it
if isfield(gt, 'elem2dof'), e2d = double(gt.elem2dof); return; end
T = double(gt.T);  edge = double(gt.F(1:2,:)).';
[~, e1] = ismember(sort(T([2 3],:).',2), edge, 'rows');
[~, e2] = ismember(sort(T([1 3],:).',2), edge, 'rows');
[~, e3] = ismember(sort(T([1 2],:).',2), edge, 'rows');
e2d = [e1 e2 e3];
end
