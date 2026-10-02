function FV = fv_setup(mesh, model, faceTag)
%FV_SETUP  Mesh-only data of the finite volume scheme on the dual mesh.
%
%   FV = fv_setup(mesh, model)            all boundary faces are walls
%   FV = fv_setup(mesh, model, faceTag)   faceTag(j) = 0 interior | 1 inflow |
%                                         2 outflow | 3 wall  (dual faces)

g = mesh.g;  gt = mesh.gt;
NC = numel(g.A);  nf = size(g.F,2);
xC = g.xC(:);  yC = g.yC(:);
K1 = g.F(5,:).';  K2 = g.F(6,:).';
if nargin < 3 || isempty(faceTag)
    faceTag = 3*(K2 <= 0);
end

%% MUSCL 
rows = zeros(12*NC,1);  cols = rows;  vx = rows;  vy = rows;  nz = 0;
for i = 1:NC
    nb = g.Nb{i};  m = numel(nb);
    if m == 0, continue; end
    [~, o] = sort(atan2(yC(nb)-yC(i), xC(nb)-xC(i)));  nb = nb(o);
    idx = zeros(1,3*m);  cx = idx;  cy = idx;  accA = 0;
    for q = 1:m
        j = nb(q);  l = nb(1 + mod(q,m));                   % cyclic
        x1 = xC(i); y1 = yC(i); x2 = xC(j); y2 = yC(j); x3 = xC(l); y3 = yC(l);
        accA = accA + 0.5*abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
        r = 3*(q-1) + (1:3);
        idx(r) = [i, j, l];
        cx(r)  = 0.5*[y2-y3, y3-y1, y1-y2];
        cy(r)  = 0.5*[x3-x2, x1-x3, x2-x1];
    end
    if accA == 0, continue; end                             
    n = 3*m;
    rows(nz+(1:n)) = i;  cols(nz+(1:n)) = idx;
    vx(nz+(1:n)) = cx/accA;  vy(nz+(1:n)) = cy/accA;  nz = nz + n;
end
FV.Gx = sparse(rows(1:nz), cols(1:nz), vx(1:nz), NC, NC);
FV.Gy = sparse(rows(1:nz), cols(1:nz), vy(1:nz), NC, NC);

%% Gauss points on the dual faces (q = 1,2 stacked)
sg = [-1 1]/sqrt(3);
f  = [(1:nf).'; (1:nf).'];
th = [-g.N(2,:); g.N(1,:)];                               % unit tangent
P  = [g.M + 0.5*sg(1)*g.L.*th, g.M + 0.5*sg(2)*g.L.*th];
FV.xq = P(1,:).';  FV.yq = P(2,:).';
FV.K1 = K1(f);  FV.K2 = K2(f);
FV.nx = g.N(1,f).';  FV.ny = g.N(2,f).';  FV.ln = g.L(f).';
FV.kn = model.k(1)*FV.nx + model.k(2)*FV.ny;
FV.tag = faceTag(f);  FV.tag = FV.tag(:);
t  = mesh.itri(f);
X1 = gt.V(:, gt.T(1,t));  X2 = gt.V(:, gt.T(2,t));  X3 = gt.V(:, gt.T(3,t));
den = (X2(2,:)-X3(2,:)).*(X1(1,:)-X3(1,:)) + (X3(1,:)-X2(1,:)).*(X1(2,:)-X3(2,:));
l1 = ((X2(2,:)-X3(2,:)).*(P(1,:)-X3(1,:)) + (X3(1,:)-X2(1,:)).*(P(2,:)-X3(2,:)))./den;
l2 = ((X3(2,:)-X1(2,:)).*(P(1,:)-X3(1,:)) + (X1(1,:)-X3(1,:)).*(P(2,:)-X3(2,:)))./den;
FV.shp = 1 - 2*[l1; l2; 1-l1-l2].';                       % CR basis at the Gauss points
FV.dof = mesh.elem2dof(t,:);
FV.nf = nf;  FV.NC = NC;  FV.A = g.A(:);  FV.xC = xC;  FV.yC = yC;

%% slope limiter (3.3) and scaling limiters: stencils and control points
in = K2 > 0;
FV.limCell = [K1; K2(in)];                                % (cell, face) pairs
FV.limNb   = [K2; K1(in)];                                % cell across the face (-1 on dOmega)
FV.limDx   = g.M(1,[1:nf, find(in).']).' - xC(FV.limCell);
FV.limDy   = g.M(2,[1:nf, find(in).']).' - yC(FV.limCell);
FV.nbA = [K1(in); K2(in)];   FV.nbB = [K2(in); K1(in)];
VX = nan(NC,4);  VY = VX;                                 % control points: vertices of K
for c = 1:4
    v = g.T(c,:).';  ok = v > 0;
    VX(ok,c) = g.V(1,v(ok));  VY(ok,c) = g.V(2,v(ok));
end
FV.VX = VX;  FV.VY = VY;

%% CFL
dia  = g.T(4,:).' > 0;
m1 = FV.A;  m2 = FV.A;                                    % boundary cells: T_K1 = T_K2 = K
tri = @(a,b,c) 0.5*abs((b(1,:)-a(1,:)).*(c(2,:)-a(2,:)) - (c(1,:)-a(1,:)).*(b(2,:)-a(2,:))).';
Td = g.T(:,dia);
m1(dia) = tri(g.V(:,Td(1,:)), g.V(:,Td(2,:)), g.V(:,Td(3,:)));
m2(dia) = tri(g.V(:,Td(1,:)), g.V(:,Td(3,:)), g.V(:,Td(4,:)));
Lf = g.L(:);
perim = accumarray([K1; K2(in)], [Lf; Lf(in)], [NC 1]);
FV.dt_geo = min(min(m1, m2)./perim)/3;
end
