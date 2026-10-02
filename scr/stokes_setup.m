function S = stokes_setup(mesh, bc)
%STOKES_SETUP  Mesh-only part of the Crouzeix-Raviart/P0 Stokes solver 
%
%   S = stokes_setup(mesh)        u = 0 on the whole boundary (Examples 1-4)
%   S = stokes_setup(mesh, bc)    non-homogeneous Dirichlet  (Example 5)

node = mesh.gt.V.';  elem = mesh.gt.T.';  elem2dof = mesh.elem2dof;  edge = mesh.edge;
NE = size(edge,1);  NT = size(elem,1);  Nu = NE;  Np = NT;

% gradients of barycentric coordinates and areas
v1 = node(elem(:,3),:) - node(elem(:,2),:);
v2 = node(elem(:,1),:) - node(elem(:,3),:);
v3 = node(elem(:,2),:) - node(elem(:,1),:);
sa = 0.5*(v2(:,1).*v3(:,2) - v2(:,2).*v3(:,1));          % signed area
Dl = zeros(NT,2,3);
Dl(:,:,1) = [-v1(:,2), v1(:,1)]./(2*sa);
Dl(:,:,2) = [-v2(:,2), v2(:,1)]./(2*sa);
Dl(:,:,3) = [-v3(:,2), v3(:,1)]./(2*sa);
area = abs(sa);
G = -2*Dl;                                               % grad psi_i (CR basis)

%% viscous term  int_T mu eps(u):eps(v)  
nb = 9*NT;  ii = zeros(nb,1);  jj = ii;  b11 = ii;  b22 = ii;  b12 = ii;  b21 = ii;  r0 = 0;
for i = 1:3
    for j = 1:3
        r = r0 + (1:NT);  r0 = r0 + NT;
        ii(r) = elem2dof(:,i);  jj(r) = elem2dof(:,j);
        b11(r) = area.*(G(:,1,i).*G(:,1,j) + 0.5*G(:,2,i).*G(:,2,j));
        b22(r) = area.*(G(:,2,i).*G(:,2,j) + 0.5*G(:,1,i).*G(:,1,j));
        b12(r) = area.*(0.5*G(:,2,i).*G(:,1,j));
        b21(r) = area.*(0.5*G(:,1,i).*G(:,2,j));
    end
end
S.rowA = [ii; ii+Nu; ii; ii+Nu];
S.colA = [jj; jj+Nu; jj+Nu; jj];
S.valA = [b11; b22; b12; b21];
S.muA  = repmat((1:NT).', 36, 1);

%% jump stabilisation  sum_e mu_e/h_e int_e [u].[v] ds  
tt = repmat((1:NT).', 3, 1);
[es, o] = sort(elem2dof(:));  ts = tt(o);
first = [true; diff(es) ~= 0];
T1 = ts(first);  T2 = T1;
isInt = accumarray(es, 1, [NE 1]) == 2;
fi = find(first);  T2(isInt) = ts(fi(isInt) + 1);
Pa = node(edge(:,1),:);  Pb = node(edge(:,2),:);
he = sqrt(sum((Pb-Pa).^2, 2));
gp = [0.5-0.5/sqrt(3); 0.5+0.5/sqrt(3)];
ldof = [elem2dof(T1,:), elem2dof(T2,:)];
tr = zeros(NE,6,2);
for q = 1:2
    Xq = Pa + gp(q)*(Pb - Pa);
    for i = 1:3
        tr(:,i,q)   =   1 - 2*(1 + sum((Xq - node(elem(T1,i),:)).*Dl(T1,:,i), 2));
        tr(:,3+i,q) = -(1 - 2*(1 + sum((Xq - node(elem(T2,i),:)).*Dl(T2,:,i), 2)));
    end
end
tr(~isInt,4:6,:) = 0;                                  % on dOmega the jump is the trace
nj = 36*NE;  iJ = zeros(nj,1);  jJ = iJ;  bJ = iJ;  r0 = 0;
for a = 1:6
    for b = 1:6
        r = r0 + (1:NE);  r0 = r0 + NE;
        iJ(r) = ldof(:,a);  jJ(r) = ldof(:,b);
        bJ(r) = 0.5*(tr(:,a,1).*tr(:,b,1) + tr(:,a,2).*tr(:,b,2));   
    end
end
S.rowJ = [iJ; iJ+Nu];  S.colJ = [jJ; jJ+Nu];  S.valJ = [bJ; bJ];
S.muJ  = repmat((1:NE).', 72, 1);
S.T1 = T1;  S.T2 = T2;

%% divergence 
rowB = repmat((1:Np).', 3, 1);
gx = G(:,1,:).*area;  gy = G(:,2,:).*area;
S.B = -[sparse(rowB, elem2dof(:), gx(:), Np, Nu), sparse(rowB, elem2dof(:), gy(:), Np, Nu)];

%% dofs, quadrature weights, boundary data
S.free = find(accumarray(elem2dof(:), 1, [Nu 1]) == 2);    % interior edges
S.free = [S.free; Nu + S.free];
S.wE   = accumarray(elem2dof(:), repmat(area/3, 3, 1), [Nu 1]);   % int_Omega psi_e = |K|
S.area = area;  S.elem2dof = elem2dof;  S.Nu = Nu;  S.Np = Np;
S.uD = zeros(2*Nu,1);  S.nit = [];
if nargin > 1 && ~isempty(bc)
    S.uD = bc.uD;  S.nit = bc.nit;
end
end
