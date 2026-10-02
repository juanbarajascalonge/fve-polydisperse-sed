function mesh = load_mesh(file)
%LOAD_MESH  Read a mesh .mat file and return the structure of BUILD_MESH.

S = load(file);
if ~(isfield(S,'g') && isfield(S,'gt'))
    mesh = build_mesh(S.Vtri, S.Ttri);
else
    g = S.g;  gt = S.gt;
    if isfield(S,'V'),    g.V  = S.V;    end
    if isfield(S,'T'),    g.T  = S.T;    end
    if isfield(S,'Vtri'), gt.V = S.Vtri; elseif isfield(gt,'Vtri'), gt.V = gt.Vtri; end
    if isfield(S,'Ttri'), gt.T = S.Ttri; elseif isfield(gt,'Ttri'), gt.T = gt.Ttri; end
    for f = {'V','T','F','A','N','M','L','xC','yC'}
        g.(f{1}) = double(g.(f{1}));
    end
    for f = {'V','T','F','A','N','M','L','B'}
        gt.(f{1}) = double(gt.(f{1}));
    end
    g.Nb        = cellfun(@(x) double(x(:).'), g.Nb,        'UniformOutput', false);
    g.cellFaces = cellfun(@(x) double(x(:).'), g.cellFaces, 'UniformOutput', false);
    gt = rmfield_if(gt, {'Vtri','Ttri'});

    % CR dofs: local edge i of a triangle is opposite vertex i (edges = gt.F(1:2,:))
    edge = gt.F(1:2,:).';
    Tt   = gt.T;
    [~, e1] = ismember(sort(Tt([2 3],:).',2), edge, 'rows');
    [~, e2] = ismember(sort(Tt([1 3],:).',2), edge, 'rows');
    [~, e3] = ismember(sort(Tt([1 2],:).',2), edge, 'rows');
    elem2dof = [e1 e2 e3];

    % primal triangle containing each dual face
    nv   = size(gt.V, 2);
    itri = max(g.F(1:2,:), [], 1).' - nv;
    bf   = g.F(6,:) <= 0;
    itri(bf) = gt.F(5, g.F(5,bf));
    if any(itri < 1 | itri > size(gt.T,2))
        error('load_mesh:itri', 'Unexpected dual-face numbering in %s.', file);
    end

    mesh.g = g;  mesh.gt = gt;  mesh.itri = itri;
    mesh.elem2dof = elem2dof;   mesh.edge = edge;
end
fprintf('mesh %s: %d triangles, %d control volumes\n', file, size(mesh.gt.T,2), numel(mesh.g.A));
end

function s = rmfield_if(s, names)
for k = 1:numel(names)
    if isfield(s, names{k}), s = rmfield(s, names{k}); end
end
end
