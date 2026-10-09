function [rowForColumn,totalCost] = assignment_min_cost(costMatrix)
%ASSIGNMENT_MIN_COST Minimum-cost assignment for a square real matrix.
% Base-MATLAB primal-dual shortest-augmenting-path implementation.
%
% rowForColumn(j) is the row assigned to column j. The implementation was
% written for this reproducibility package and requires no toolbox.

[nRows,nCols] = size(costMatrix);
if nRows ~= nCols
    error('assignment_min_cost:SquareMatrixRequired', ...
        'The cost matrix must be square.');
end
if ~isreal(costMatrix) || any(~isfinite(costMatrix(:)))
    error('assignment_min_cost:FiniteRealMatrixRequired', ...
        'All costs must be finite and real.');
end

n = nRows;
if n == 0
    rowForColumn = zeros(1,0);
    totalCost = 0;
    return;
end

% Index 1 is a dummy column. Real columns are stored at indices 2:n+1.
u = zeros(n,1);
v = zeros(n+1,1);
p = zeros(n+1,1);
way = zeros(n+1,1);

for newRow = 1:n
    p(1) = newRow;
    dummyOrColumn = 1;
    minReduced = inf(n+1,1);
    used = false(n+1,1);

    while true
        used(dummyOrColumn) = true;
        activeRow = p(dummyOrColumn);
        delta = inf;
        nextColumn = 0;

        for columnIndex = 2:(n+1)
            if ~used(columnIndex)
                reduced = costMatrix(activeRow,columnIndex-1) ...
                    - u(activeRow) - v(columnIndex);
                if reduced < minReduced(columnIndex)
                    minReduced(columnIndex) = reduced;
                    way(columnIndex) = dummyOrColumn;
                end
                if minReduced(columnIndex) < delta
                    delta = minReduced(columnIndex);
                    nextColumn = columnIndex;
                end
            end
        end

        if ~isfinite(delta) || nextColumn == 0
            error('assignment_min_cost:NoAugmentingPath', ...
                'No finite augmenting path was found.');
        end

        for columnIndex = 1:(n+1)
            if used(columnIndex)
                if p(columnIndex) ~= 0
                    u(p(columnIndex)) = u(p(columnIndex)) + delta;
                end
                v(columnIndex) = v(columnIndex) - delta;
            else
                minReduced(columnIndex) = minReduced(columnIndex) - delta;
            end
        end

        dummyOrColumn = nextColumn;
        if p(dummyOrColumn) == 0
            break;
        end
    end

    while true
        previous = way(dummyOrColumn);
        p(dummyOrColumn) = p(previous);
        dummyOrColumn = previous;
        if dummyOrColumn == 1
            break;
        end
    end
end

rowForColumn = p(2:end).';
linearIndex = sub2ind([n n],rowForColumn,1:n);
totalCost = sum(costMatrix(linearIndex));
end
