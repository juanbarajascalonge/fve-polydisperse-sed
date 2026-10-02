function ex = manufactured_ex1(model)
%MANUFACTURED_EX1  Exact solution and forcing terms of Example 1 
%
%   ex = manufactured_ex1(model)   returns function handles
%     ex.phi(x,y,t)    N columns,  Phi_ex = phi0 * (1/2, 1/3, 1/6)
%     ex.u(x,y)        [u1 u2],    ex.gradu(x,y) = [u1_x u1_y u2_x u2_y]
%     ex.p(x,y)
%     ex.source(x,y,t) forcing of the transport equations (N columns)
%     ex.force(x,y,t)  extra force of the Stokes problem
%                      -div(mu(phi_ex) eps(u_ex)) + grad p_ex - g(phi_ex)
c = [1/2 1/3 1/6];
ex.phi    = @(x,y,t) phi0(x,y,t)*c;
ex.u      = @vel;
ex.gradu  = @gradvel;
ex.p      = @(x,y) cos(pi*x) + cos(pi*y);
ex.source = @(x,y,t) source(x,y,t,c,model);
ex.force  = @(x,y,t) force(x,y,t,model);
end

function s = phi0(x,y,t)
s = 0.2*sin(pi*x).^2.*sin(pi*y).^2*cos(t);
end

function [st, sx, sy] = dphi0(x,y,t)
st = -0.2*sin(pi*x).^2.*sin(pi*y).^2*sin(t);
sx =  0.2*pi*sin(2*pi*x).*sin(pi*y).^2*cos(t);
sy =  0.2*pi*sin(pi*x).^2.*sin(2*pi*y)*cos(t);
end

function U = vel(x,y)
U = [sin(2*pi*y).*(1 - cos(2*pi*x)), sin(2*pi*x).*(cos(2*pi*y) - 1)];
end

function D = gradvel(x,y)
D = 2*pi*[ sin(2*pi*x).*sin(2*pi*y), cos(2*pi*y).*(1 - cos(2*pi*x)), ...
           cos(2*pi*x).*(cos(2*pi*y) - 1), -sin(2*pi*x).*sin(2*pi*y)];
end

function S = source(x,y,t,c,model)
% S_l = d_t phi_l + u.grad phi_l + k.grad f_l(Phi_ex)        (div u = 0)
s = phi0(x,y,t);  [st, sx, sy] = dphi0(x,y,t);
U = vel(x,y);
n = model.nrz;  d = model.delta;
Vl  = (1 - s).^(n-2);                  % M x N
dVl = -(n-2).*(1 - s).^(n-3);
G   = sum(d.*c.*Vl, 2);   dG = sum(d.*c.*dVl, 2);
v   = (1 - s).*(d.*Vl - s.*G);                                  % v_l(c*s)
dv  = -(d.*Vl - s.*G) + (1 - s).*(d.*dVl - G - s.*dG);          % d v_l / ds
dF  = c.*(v + s.*dv);                                           % d f_l / ds
S   = c.*(st + U(:,1).*sx + U(:,2).*sy) + dF.*(model.k(1)*sx + model.k(2)*sy);
end

function F = force(x,y,t,model)
s = phi0(x,y,t);  [~, sx, sy] = dphi0(x,y,t);
D = gradvel(x,y);                      % u1x u1y u2x u2y
lap1 = 4*pi^2*sin(2*pi*y).*(2*cos(2*pi*x) - 1);
lap2 = 4*pi^2*sin(2*pi*x).*(1 - 2*cos(2*pi*y));
mu   = (1/model.mu0)*(1 - s/model.phimax).^(-model.upsilon);
dmu  = (1/model.mu0)*model.upsilon*(1 - s/model.phimax).^(-model.upsilon-1)/model.phimax;
mux  = dmu.*sx;  muy = dmu.*sy;
exy  = 0.5*(D(:,2) + D(:,3));
% -div(mu eps(u)) = -(mu/2) Lap u - eps(u) grad mu
F = [-0.5*mu.*lap1 - (D(:,1).*mux + exy.*muy), ...
     -0.5*mu.*lap2 - (exy.*mux + D(:,4).*muy)];
F = F + [-pi*sin(pi*x), -pi*sin(pi*y)] - model.dr*s.*model.k;
end
