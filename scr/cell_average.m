function U = cell_average(f, mesh)
%CELL_AVERAGE  Averages of f over the dual control volumes.
%
%   U = cell_average(f, mesh)
%   f : @(x,y) returning an M x N array for column vectors x, y (M points)
%   U : NC x N

g = mesh.g;  T = g.T;  V = g.V;  NC = size(T,2);
gl = [0.069431844202974 0.330009478207572 0.669990521792428 0.930568155797026];
wl = [0.173927422568727 0.326072577431273 0.326072577431273 0.173927422568727];
[s, r] = ndgrid(gl, gl);   [ws, wr] = ndgrid(wl, wl);
l1 = s(:);  l2 = r(:).*(1 - s(:));  w = ws(:).*wr(:).*(1 - s(:))*2;   % sum(w) = 1

U = 0;
for part = 1:2
    if part == 1
        c = 1:NC;            P = T([1 2 3], :);
    else
        c = find(T(4,:) > 0); P = T([1 3 4], c);
    end
    A = V(:,P(1,:));  B = V(:,P(2,:));  C = V(:,P(3,:));
    ar = 0.5*abs((B(1,:)-A(1,:)).*(C(2,:)-A(2,:)) - (C(1,:)-A(1,:)).*(B(2,:)-A(2,:)));
    I = 0;
    for q = 1:numel(w)
        X = l1(q)*A + l2(q)*B + (1 - l1(q) - l2(q))*C;
        I = I + w(q)*f(X(1,:).', X(2,:).');
    end
    if part == 1, U = zeros(NC, size(I,2)); end
    U(c,:) = U(c,:) + ar(:).*I;
end
U = U./g.A(:);
end
