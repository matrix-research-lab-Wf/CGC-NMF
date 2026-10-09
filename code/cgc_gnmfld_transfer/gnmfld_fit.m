function V = gnmfld_fit(X, y, labeled_idx, c, W, alpha, beta, U0, V0, max_iter, eps0)
%GNMFLD_FIT Memory-efficient multiplicative solver for GNMFLD.
% MATLAB R2019a compatible.
%
% Objective:
%   ||X-UV'||_F^2
%   + alpha ||P_Omega(V)-P_Omega(Y)||_F^2
%   + beta Tr(V' L V)
%
% Memory strategy:
% - do NOT form a full n-by-c dense label matrix Y;
% - do NOT form full XtU, numerator, and denominator simultaneously;
% - compute W*V once per V-update, then update V in row blocks;
% - the row-block update is algebraically identical to the simultaneous
%   multiplicative update because W*V is computed from the old V before
%   any V row is overwritten.

n = size(X,2);
m = size(X,1);

U = U0;
V = V0;

if size(U,1) ~= m || size(U,2) ~= c
    error('U0 has incompatible dimensions.');
end
if size(V,1) ~= n || size(V,2) ~= c
    error('V0 has incompatible dimensions.');
end

% label_class(i)=k for labeled sample i in class k, else 0.
label_class = zeros(n,1);
label_class(labeled_idx) = y(labeled_idx);

D = full(sum(W,2));

% Smaller blocks reduce peak memory on COIL100/MNIST.
block_size = 256;

for it = 1:max_iter
    % ---- U update ----
    VtV = V' * V;                  % c-by-c
    numU = X * V;                  % m-by-c
    denU = U * VtV;                % m-by-c
    denU(denU < eps0) = eps0;
    U = U .* (numU ./ denU);

    clear VtV numU denU

    % ---- V update ----
    % Compute graph contribution from OLD V once. This lets us overwrite
    % V block by block without changing the intended simultaneous update.
    WV = W * V;                    % n-by-c
    UtU = U' * U;                  % c-by-c

    for first = 1:block_size:n
        last = min(n, first + block_size - 1);
        rows = first:last;
        nb = numel(rows);

        Vb = V(rows,:);

        % Reconstruction numerator for this block only.
        num = X(:,rows)' * U;      % nb-by-c
        num = num + beta * WV(rows,:);

        % Add alpha*Y only at labeled rows; Y is one-hot.
        lbl = label_class(rows);
        pos = find(lbl > 0);
        if ~isempty(pos)
            lin = sub2ind([nb,c], pos, lbl(pos));
            num(lin) = num(lin) + alpha;
        end

        % Denominator for this block only.
        den = Vb * UtU;
        den = den + beta * bsxfun(@times, D(rows), Vb);

        % Add alpha*P_Omega(V) only on labeled rows.
        if ~isempty(pos)
            den(pos,:) = den(pos,:) + alpha * Vb(pos,:);
        end

        den(den < eps0) = eps0;
        V(rows,:) = Vb .* (num ./ den);

        clear Vb num den lbl pos lin
    end

    clear WV UtU
end
end
