function [Acomplete,difference] = AltMin3(AA,MASK,r,maxIter,tol,show)

% This is an modified Alternative Minimization from AltMin2, replacing CVX
% with lsqr, with the purpose of reducing time cost.

% Iteration log will show if show = true.

[m,n] = size(AA);
NumOfSample = sum(MASK(:));

ROW = Row_Nonzero(M); % Nonzeros entry index for each row
COLUMN = Column_Nonzero(M); % Nonzeros entry index for each column

[U,~,~] = svd(M);
U = U(:,1:r);

%U = OrthoNormal(U);

difference = [];

for t = 1:maxIter
    if show
        display(t);
    end 
    V = zeros(r,n);
    for j = 1:n
        pos = cut_zero(COLUMN(:,j));
        V(:,j) = U(pos,:) \ AA(pos,j);
%         V(:,j) = pinv(U(pos,:)) * AA(pos,j);
%         V(:,j) = lsqr(U(pos,:),AA(pos,j));    
        if r == 1 && isinf(V(j))
            V(j) = 0;
        end
    end
%     temp = norm((U*V-A),'fro')/norm(A,'fro');
    temp = norm((MASK.*(U*V-AA)),'fro')/norm(AA,'fro'); % difference on masked samples
    difference = [difference temp];
    if temp < tol
        break;
    end
    
    U = zeros(m,r);
    for i = 1:m
        pos = cut_zero(ROW(i,:));
        u = V(:,pos)' \ AA(i,pos)';
%         u = pinv(V(:,pos)') * AA(i,pos)';
%         u = lsqr(V(:,pos)',AA(i,pos)');
        U(i,:) = u';
        if r == 1 && isinf(U(i))
            U(i) = 0;
        end
    end
%     temp = norm((U*V-A),'fro')/norm(A,'fro');
    temp = norm((MASK.*(U*V-AA)),'fro')/norm(AA,'fro'); % difference on masked samples
    difference = [difference temp];
    if temp < tol
        break;
    end
    
end

Acomplete = U*V;
% Acomplete = Acomplete .* ~MASK;
% Acomplete = Acomplete + AA;