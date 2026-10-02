function [R, alpha] = fv_residual(Phi, ux, uy, FV, model, irp, Phi_in)
%FV_RESIDUAL  Right-hand side of the semi-discrete FV scheme (3.4).
%
%   [R, alpha] = fv_residual(Phi, ux, uy, FV, model, irp, Phi_in)
%
%   Phi    : NC x N cell averages,   ux,uy : CR velocity dofs
%   irp    : true -> MUSCL + scaling limiters (IRP-FVE), false -> MUSCL with the
%            slope limiter only (FVE); 
%   Phi_in : 1 x N feed concentration on inflow faces (tag 1), optional
%   R      : dPhi/dt = R,   alpha : max_{K,L} alpha_KL 

if nargin < 7, Phi_in = []; end
N = size(Phi,2);  nf = FV.nf;

% piecewise linear reconstruction  P_K(x) = a0 + a1*x + a2*y
[a0, a1, a2] = reconstruct(Phi, FV, model.phimax, irp);

% states on both sides of every Gauss point
uL = a0(FV.K1,:) + a1(FV.K1,:).*FV.xq + a2(FV.K1,:).*FV.yq;
inner = FV.tag == 0;
uR = uL;
uR(inner,:) = a0(FV.K2(inner),:) + a1(FV.K2(inner),:).*FV.xq(inner) + a2(FV.K2(inner),:).*FV.yq(inner);
inflow = FV.tag == 1;  outflow = FV.tag == 2;  wall = FV.tag == 3;
if any(inflow), uR(inflow,:) = repmat(Phi_in(:).', nnz(inflow), 1); end

% normal velocity u_h . n from the CR basis of the containing triangle
qn = sum(ux(FV.dof).*FV.shp, 2).*FV.nx + sum(uy(FV.dof).*FV.shp, 2).*FV.ny;

% LLF flux 
[fn, al] = llf(uL, uR, qn, FV.kn, model);
fn(wall,:)    = 0;                                   % zero total flux (1.3)
fn(inflow,:)  = qn(inflow).*uR(inflow,:);            % upwind: feed state
fn(outflow,:) = qn(outflow).*uL(outflow,:);          % upwind: interior state
al(inflow | outflow) = abs(qn(inflow | outflow));
alpha = max(al(~wall));

% assemble (two-point Gauss rule on each face)
c  = 0.5*FV.ln(1:nf).*(fn(1:nf,:) + fn(nf+1:end,:));
K1 = FV.K1(1:nf);  K2 = FV.K2(1:nf);  in = K2 > 0;
R = zeros(FV.NC, N);
for l = 1:N
    R(:,l) = -accumarray(K1, c(:,l), [FV.NC 1]) + accumarray(K2(in), c(in,l), [FV.NC 1]);
end
R = R./FV.A;
end

%% ========================================================================
function [a0, a1, a2] = reconstruct(Phi, FV, phimax, irp)
h = 1;  b = 1;  eps_pos = 1e-100;

[NC, N] = size(Phi);  xC = FV.xC;  yC = FV.yC;
a1 = FV.Gx*Phi;  a2 = FV.Gy*Phi;                 

%% (1) MUSCL slope limiter
pmin = Phi;  pmax = Phi;
for l = 1:N
    pmin(:,l) = min(Phi(:,l), accumarray(FV.nbA, Phi(FV.nbB,l), [NC 1], @min,  inf));
    pmax(:,l) = max(Phi(:,l), accumarray(FV.nbA, Phi(FV.nbB,l), [NC 1], @max, -inf));
end
c  = FV.limCell;  j = FV.limNb;                  % (cell, face) pairs and neighbour
du = a1(c,:).*FV.limDx + a2(c,:).*FV.limDy;      % Phi*_{K,L} - Phi_K
r  = ones(size(du));
up = du > 1e-14;  dn = du < -1e-14;
dmax = pmax(c,:) - Phi(c,:);  dmin = pmin(c,:) - Phi(c,:);
r(up) = dmax(up)./du(up);
r(dn) = dmin(dn)./du(dn);
uij = max(0, min(min(1, b*r), b));
aij = uij;                                       % boundary faces
in  = j > 0;
if h < 1
    sx  = xC(j(in)) - xC(c(in));  sy = yC(j(in)) - yC(c(in));
    den = Phi(j(in),:) - Phi(c(in),:);
    r1  = ones(size(den));  r2 = r1;  ok = abs(den) >= 1e-14;
    gi  = a1(c(in),:).*sx + a2(c(in),:).*sy;     % u_i - W_i*
    gj  = a1(j(in),:).*sx + a2(j(in),:).*sy;     % W_j* - u_j
    r1(ok) = gi(ok)./den(ok);  r2(ok) = gj(ok)./den(ok);
    aij(in,:) = min(max(0, min(1, r1)), max(0, min(1, r2)));
end
tij = h*uij + (1-h)*aij;
th  = ones(NC, N);
for l = 1:N
    th(:,l) = accumarray(c, tij(:,l), [NC 1], @min, 1);
end
a1 = a1.*th;  a2 = a2.*th;

if irp
    %% (2A) positivity of each component at the vertices of K
    for l = 1:N
        m   = min(Phi(:,l) + a1(:,l).*(FV.VX - xC) + a2(:,l).*(FV.VY - yC), [], 2);
        den = Phi(:,l) - m;
        t   = ones(NC,1);
        bad = m < eps_pos;
        t(bad & den > 0)  = min(1, (Phi(bad & den > 0,l) - eps_pos)./den(bad & den > 0));
        t(bad & den <= 0) = 0;
        a1(:,l) = a1(:,l).*t;  a2(:,l) = a2(:,l).*t;
    end
    %% (2B) upper bound of the sum at the vertices of K
    phi = sum(Phi, 2);
    M   = max(phi + sum(a1,2).*(FV.VX - xC) + sum(a2,2).*(FV.VY - yC), [], 2);
    den = M - phi;
    t   = ones(NC,1);
    bad = M > phimax;
    t(bad & den > 0)  = max(0, min(1, (phimax - phi(bad & den > 0))./den(bad & den > 0)));
    t(bad & den <= 0) = 0;
    a1 = a1.*t;  a2 = a2.*t;
end
a0 = Phi - a1.*xC - a2.*yC;
end

%% ========================================================================
function [f, alpha] = llf(uL, uR, qn, kn, model)
% local Lax-Friedrichs flux for F(U).n = (u.n) U + (k.n) f(U); the wave
% speeds are bounded by (2.9) with the bounds M1 <= lambda, v <= M2 of (1.7)
[fL, M1L, M2L] = mlb(uL, model);
[fR, M1R, M2R] = mlb(uR, model);
alpha = max(abs(qn + kn.*min(M1L, M1R)), abs(qn + kn.*max(M2L, M2R)));
f = 0.5*(qn.*(uL + uR) + kn.*(fL + fR) - alpha.*(uR - uL));
end

%% ========================================================================
function [f, M1, M2] = mlb(U, model)
% MLB flux
n = model.nrz;  d = model.delta;  pm = model.phimax;
phi = sum(U, 2);
ok  = phi >= 0 & phi <= pm;
s   = min(max(phi, 0), pm);
ps  = ((n-2)*pm - 1)./(n-3);                       % phi_{l,*}
V   = (1-ps).^(n-2) - (n-2).*(1-ps).^(n-3).*(s - ps);     % tangent branch
dV  = -(n-2).*(1-ps).^(n-3).*ones(size(s));
pw  = s < ps;                                      % power branch
Vp  = (1-s).^(n-2);   dVp = -(n-2).*(1-s).^(n-3);
V(pw) = Vp(pw);  dV(pw) = dVp(pw);
V(~ok,:) = 0;  dV(~ok,:) = 0;

W  = d.*V;
B  = W - sum(U.*W, 2);                             
f  = (1 - phi).*U.*B;
M1 = sum(d.*U.*((1-phi).*dV - 2*V), 2);
M2 = (1 - phi).*B(:,1);
end
