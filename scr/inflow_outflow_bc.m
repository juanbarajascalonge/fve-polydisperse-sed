function bc = inflow_outflow_bc(mesh, geom, Ubar, Phi_in)
%INFLOW_OUTFLOW_BC  Boundary data of Example 5 (Section 4.5).
%
%   bc = inflow_outflow_bc(mesh, geom, Ubar, Phi_in)
%
%   geom.x0, geom.x1     abscissae of the walls holding Gamma_in and Gamma_out
%   geom.yin, geom.yout  [y1 y2] of the two (vertical) mouths;  geom.tol
%   Ubar                 mean inflow velocity
%   Phi_in               1 x N feed concentration
%
%   bc.faceTag(j)  0 interior | 1 inflow | 2 outflow | 3 wall   (dual faces)
%   bc.edgeTag(e)  idem for primal edges (= dual cells = CR dofs)
%   bc.uD          CR Dirichlet data: 
%   bc.nit         u_D part of the jump term on the mouths: with u = u_D ~= 0
%                  the term int_e [u_h].[v_h] must read int_e (u_h - u_D).v_h

g = mesh.g;  gt = mesh.gt;  tol = geom.tol;
nf = size(g.F,2);  NE = size(gt.F,2);

%% tags
bc.faceTag = zeros(nf,1);  bc.edgeTag = zeros(NE,1);
for j = find(g.F(6,:) <= 0)
    xm = g.M(1,j);  ym = g.M(2,j);
    if abs(xm - geom.x0) < tol && ym > geom.yin(1) - tol && ym < geom.yin(2) + tol
        tg = 1;
    elseif abs(xm - geom.x1) < tol && ym > geom.yout(1) - tol && ym < geom.yout(2) + tol
        tg = 2;
    else
        tg = 3;
    end
    bc.faceTag(j) = tg;  bc.edgeTag(g.F(5,j)) = tg;      % dual cell == primal edge
end
Lin  = sum(gt.L(bc.edgeTag == 1));  Lout = sum(gt.L(bc.edgeTag == 2));
if abs(Lin - diff(geom.yin)) > 1e-8 || abs(Lout - diff(geom.yout)) > 1e-8
    error('inflow_outflow_bc:mouths', 'The end points of the mouths must be mesh nodes.');
end

%% Dirichlet data and Nitsche weights
gp = [0.5-0.5/sqrt(3); 0.5+0.5/sqrt(3)];  prof = @(s) 6*s.*(1-s);
E  = find(bc.edgeTag == 1 | bc.edgeTag == 2);  nb = numel(E);
uD = zeros(2*NE,1);
nit.tri = gt.F(5,E).';  nit.ldof = mesh.elem2dof(nit.tri,:);
nit.Wx = zeros(nb,3);  nit.Wy = zeros(nb,3);
for m = 1:nb
    e = E(m);
    a = gt.V(:,gt.F(1,e));  b = gt.V(:,gt.F(2,e));  n = gt.N(:,e);   % outward normal
    if bc.edgeTag(e) == 1, seg = geom.yin;  d = -n;  else, seg = geom.yout;  d = n;  end
    X = [gt.V(:, gt.T(:,nit.tri(m))); 1 1 1];
    for q = 1:2
        P   = a + gp(q)*(b - a);
        uq  = Ubar*prof((P(2) - seg(1))/diff(seg))*d;
        tr  = 1 - 2*(X \ [P; 1]).';                       % CR basis on the edge
        uD(e) = uD(e) + 0.5*uq(1);   uD(NE+e) = uD(NE+e) + 0.5*uq(2);
        nit.Wx(m,:) = nit.Wx(m,:) + 0.5*tr*uq(1);
        nit.Wy(m,:) = nit.Wy(m,:) + 0.5*tr*uq(2);
    end
end

%% exact discrete compatibility
Q   = @(I) sum(gt.L(I).' .* (uD(I).*gt.N(1,I).' + uD(NE+I).*gt.N(2,I).'));
in  = find(bc.edgeTag == 1);  out = find(bc.edgeTag == 2);
sc  = -Q(in)/Q(out);
uD([out; NE+out]) = sc*uD([out; NE+out]);
o = bc.edgeTag(E) == 2;
nit.Wx(o,:) = sc*nit.Wx(o,:);  nit.Wy(o,:) = sc*nit.Wy(o,:);

bc.uD = uD;  bc.nit = nit;  bc.Phi_in = Phi_in(:).';
fprintf('inflow_outflow_bc: Q_in = %.6f, Q_out = %.6f\n', -Q(in), Q(out));
end
