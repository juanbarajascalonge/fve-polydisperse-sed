function [ux, uy, p, info] = stokes_solve(S, model, Phi, F)
%STOKES_SOLVE  
%
%   -div(mu(phi) eps(u)) + grad p = g(phi) + F,   div u = 0,   int p = 0
%
%   Phi : NC x N cell averages on the dual mesh (NC = number of primal edges)
%   F   : optional NE x 2 extra force at the edge midpoints (Example 1)
%   ux,uy : CR dofs (edge means),  p : P0 pressure,
%   info.divmax = max_T |int_T div u_h|   (local DDF property)

Nu = S.Nu;  Np = S.Np;
phi  = sum(Phi, 2);
phiT = mean(phi(S.elem2dof), 2);
s    = max(model.eta, min(phiT/model.phimax, 1 - model.eta));
muT  = (1/model.mu0)*(1 - s).^(-model.upsilon);
muE  = 0.5*(muT(S.T1) + muT(S.T2));                   % average over the elements sharing e

A = sparse(S.rowA, S.colA, muT(S.muA).*S.valA, 2*Nu, 2*Nu) ...
  + sparse(S.rowJ, S.colJ, model.gamma*muE(S.muJ).*S.valJ, 2*Nu, 2*Nu);

% g(phi) = dr*phi*k; the CR mass matrix is diagonal, int psi_e psi_e = |K_e|
f = [S.wE.*(model.dr*model.k(1)*phi); S.wE.*(model.dr*model.k(2)*phi)];
if nargin > 3 && ~isempty(F)
    f = f + [S.wE.*F(:,1); S.wE.*F(:,2)];
end
if ~isempty(S.nit)            % u_D part of the jump term on the open boundaries
    sc = model.gamma*muT(S.nit.tri);
    f = f + [accumarray(S.nit.ldof(:), reshape(sc.*S.nit.Wx,[],1), [Nu 1]); ...
             accumarray(S.nit.ldof(:), reshape(sc.*S.nit.Wy,[],1), [Nu 1])];
end

% saddle point system
fr = S.free;  nf = numel(fr);  B = S.B(1:Np-1,:);
K = [A(fr,fr),  B(:,fr).';
     B(:,fr),   sparse(Np-1,Np-1)];
sol = K \ [f(fr) - A(fr,:)*S.uD; -B*S.uD];
p = [sol(nf+1:end); 0];
p = p - (S.area.'*p)/sum(S.area);
u = S.uD;  u(fr) = sol(1:nf);
ux = u(1:Nu);  uy = u(Nu+1:end);
info.divmax = max(abs(S.B*u));
end
