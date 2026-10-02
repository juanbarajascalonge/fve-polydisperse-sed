function mesh = build_mesh(Vtri, Ttri)
%BUILD_MESH  
%
%   mesh = build_mesh(Vtri, Ttri)
%     Vtri : 2 x nv vertex coordinates,  Ttri : 3 x nt triangles.
%
%   Output fields:
%     mesh.gt  primal mesh: V T F A N M L B   (F = [v1;v2;opp1;opp2;t1;t2])
%     mesh.g   dual mesh  : V T F A N M L Nb cellFaces xC yC
%              (F = [p1;p2;ref1;ref2;K1;K2], N points from K1 to K2, K2=-1 on dOmega)
%     mesh.elem2dof (nt x 3)  local edge i of a triangle is opposite vertex i
%     mesh.itri     (nf x 1)  primal triangle containing each dual face

Vtri = double(Vtri);
Ttri = sort(double(Ttri), 1);
nv = size(Vtri, 2);  nt = size(Ttri, 2);

%% ---------------- primal edges ------------------------------------------
loc = [2 3; 1 3; 1 2];                        % local edge i is opposite vertex i
E   = [Ttri(loc(1,:),:), Ttri(loc(2,:),:), Ttri(loc(3,:),:)].';   % 3nt x 2
[edge, ~, j] = unique(E, 'rows');             % lexicographic, edge(:,1) < edge(:,2)
NE  = size(edge, 1);
elem2dof = reshape(j, nt, 3);

tt  = repmat((1:nt).', 3, 1);                 % triangle of each row of E
opp = [Ttri(1,:), Ttri(2,:), Ttri(3,:)].';    % vertex opposite to each row of E
cnt = accumarray(j, 1, [NE 1]);
t1  = accumarray(j, tt, [NE 1], @min);
t2  = accumarray(j, tt, [NE 1], @max);  t2(cnt == 1) = -1;
o1  = zeros(NE,1);  o2 = -ones(NE,1);
r1  = tt == t1(j);   o1(j(r1)) = opp(r1);
r2  = tt == t2(j);   o2(j(r2)) = opp(r2);

X  = reshape(Vtri(:, Ttri), 2, 3, nt);
At = 0.5*abs( (X(1,2,:)-X(1,1,:)).*(X(2,3,:)-X(2,1,:)) ...
            - (X(1,3,:)-X(1,1,:)).*(X(2,2,:)-X(2,1,:)) );
gt.V = Vtri;   gt.T = Ttri;
gt.F = [edge.'; o1.'; o2.'; t1.'; t2.'];
gt.A = reshape(At, 1, nt);
gt.B = reshape(mean(X, 2), 2, nt);
[gt.N, gt.M, gt.L] = face_geometry(Vtri, edge(:,1), edge(:,2), o1);

%% ---------------- dual cells --------------------------------------------
V  = [Vtri, gt.B];
c2 = nv + t2;  c2(t2 < 0) = -1;
T  = [edge(:,1).'; nv + t1.'; edge(:,2).'; c2.'];
Ac = gt.A(t1);  in = t2 > 0;
Ac(in) = Ac(in) + gt.A(t2(in));
g.V = V;  g.T = T;  g.A = Ac/3;

%% ---------------- dual faces --------------------------------------------
% interior faces: [vertex v of triangle t, barycentre of t], shared by the
% two dual cells of the edges of t that contain v.
P1 = [];  P2 = [];  K1 = [];  K2 = [];  R1 = [];  R2 = [];  itri = [];
for i = 1:3
    o  = setdiff(1:3, i);
    v  = Ttri(i,:).';
    ea = elem2dof(:, o(1));   eb = elem2dof(:, o(2));
    wa = Ttri(o(2),:).';      wb = Ttri(o(1),:).';   % other endpoint of ea / eb
    sw = ea > eb;                                     % K1 = smaller cell index
    [ea(sw), eb(sw)] = deal(eb(sw), ea(sw));
    [wa(sw), wb(sw)] = deal(wb(sw), wa(sw));
    P1 = [P1; v];  P2 = [P2; nv + (1:nt).'];         
    K1 = [K1; ea]; K2 = [K2; eb];  R1 = [R1; wa];  R2 = [R2; wb];   
    itri = [itri; (1:nt).'];                          
end
% boundary faces
bE = find(~in);
P1 = [P1; edge(bE,1)];  P2 = [P2; edge(bE,2)];  R1 = [R1; nv + t1(bE)];  R2 = [R2; -ones(numel(bE),1)];
K1 = [K1; bE];          K2 = [K2; -ones(numel(bE),1)];
itri = [itri; t1(bE)];

g.F = [P1.'; P2.'; R1.'; R2.'; K1.'; K2.'];
[g.N, g.M, g.L] = face_geometry(V, P1, P2, R1);
nf = numel(K1);  NC = NE;
inF = K2 > 0;
cf  = [K1; K2(inF)];
fid = [(1:nf).'; find(inF)];
nbr = [K2(inF); K1(inF)];
g.cellFaces = accumarray(cf, fid, [NC 1], @(x) {sort(x).'}).';
g.Nb        = accumarray([K1(inF); K2(inF)], nbr, [NC 1], @(x) {unique(x).'}).';
g.xC = gt.M(1,:);   g.yC = gt.M(2,:);

mesh.g = g;  mesh.gt = gt;  mesh.itri = itri;
mesh.elem2dof = elem2dof;  mesh.edge = edge;
end

%% ------------------------------------------------------------------------
function [N, M, L] = face_geometry(V, p1, p2, ref)
a = V(:, p1);  b = V(:, p2);  d = a - b;
L = sqrt(sum(d.^2, 1));
N = [d(2,:); -d(1,:)] ./ L;
flip = sum(N .* (V(:, ref) - a), 1) > 0;
N(:, flip) = -N(:, flip);
M = 0.5*(a + b);
end
