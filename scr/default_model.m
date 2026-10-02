function model = default_model(varargin)
%DEFAULT_MODEL  Model parameters of Section 4, with name/value overrides.
%
%   model = default_model('d',[2.9 2.0]*1e-3, 'nrz',4.6, 'phimax',0.6, ...
%                         'mu0',100, 'k',[1 0])
%
%   Fields (dimensionless variables, X_ref = 1 m, T_ref = X_ref/v_St):
%     d       particle diameters d_1 >= ... >= d_N [m]   (or give 'delta')
%     delta   delta_l = d_l^2/d_1^2
%     nrz     Richardson-Zaki exponents n_l (scalar or 1 x N)
%     phimax  maximum total concentration
%     k       downward unit vector
%     mu0     viscosity parameter:  mu(phi) = (1/mu0)(1-max(eta,min(phi/phimax,1-eta)))^(-upsilon)
%     upsilon, eta, rhos, rhof  (g(phi) = (rhos-rhof)/rhof * phi * k)
%     gamma   jump-stabilisation parameter of the Stokes solver

model = struct('d',[], 'delta',[], 'nrz',4.6, 'phimax',0.6, 'k',[1 0], ...
               'mu0',1, 'upsilon',2, 'eta',1e-8, 'rhos',2790, 'rhof',1208, ...
               'gamma',1);
for i = 1:2:numel(varargin)
    if ~isfield(model, varargin{i})
        error('default_model:field', 'Unknown parameter ''%s''.', varargin{i});
    end
    model.(varargin{i}) = varargin{i+1};
end
if isempty(model.delta)
    model.delta = (model.d(:).'/model.d(1)).^2;
end
model.N     = numel(model.delta);
model.nrz   = model.nrz(:).' .* ones(1, model.N);
model.k     = model.k(:).';
model.dr    = (model.rhos - model.rhof)/model.rhof;
end
