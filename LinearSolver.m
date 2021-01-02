classdef LinearSolver
    methods
        function x = Gauss_Seidel(obj, A, b, eps, max_iter, initial_guess)
            [n,m] = size(A);
            A = [A b];
            x = initial_guess;
            k = 1;
            while k <= max_iter
                err = 0;
                for i = 1 : n
                    s = 0;
                    for j = 1: n
                        s = s - A(i,j)*x(j);
                    end
                    s = (s+A(i,n+1))/A(i,i);
                    if abs(s) > err
                        err = abs(s);
                    end
                    x(i) = x(i) + s;
                end
                if err <= eps
                    break;
                else
                    k = k+1;
                end
            end
        end
        function x = LU_Decomp(obj, A, b)
            [L,U,P] = objLU(A,0);
            y = obj.forward_sub(L, P*b);
            [n,m] = size(U);
            x = obj.back_sub(U,n,y);
        end
        
        function [L,U,P]= LU(obj, A,threshold)
            [nRow, nCol] = size(A);
            
            P=diag(ones(nRow,1));
            U=zeros(nRow);
            L=zeros(nRow);
            
            for n=1:nRow-1
                currentPivot=A(n,n);
                
                E=diag(ones(nRow,1));
                
                maxPivot=max(A(n+1:end,n));
                if abs(currentPivot)<eps  %zero, do row exchage always
                    if abs(maxPivot)<eps  % not possible to exchange
                        error 'unable to complete LU decomposition, bad A'
                    else
                        [A, E, L] = obj.flipRows(A,n);
                    end
                else %not a zero pivot, but still can exchange, check threshold
                    if abs(currentPivot)<abs(maxPivot)
                        if abs(currentPivot-maxPivot)>=threshold
                            [A, E, L] = obj.flipRows(A,n);
                        end
                    end
                end
                
                P=P*E;  %update the perumtation matrix
                
                for i=n+1:nRow
                    L(i,n)=A(i,n)/A(n,n);
                    A(i,n)=0;
                    for j=n+1:nRow
                        A(i,j)=A(i,j)-L(i,n)*A(n,j);
                    end
                end
            end
            
            L=L+diag(ones(nRow,1));
            P=P';
            U=A;
        end
        function [A,E,L] = flipRows(obj,A, n)
            [c,I]= max(abs(A(n:end,n)));
            I=I+(n-1);
            tmp=A(n,:);
            A(n,:)=A(I,:);
            A(I,:)=tmp;
            
            %make sure we also flip the L matrix rows to keep in sync
            tmp=L(n,:);
            L(n,:)=L(I,:);
            L(I,:)=tmp;
            
            %now make the elementary matrix for this move
            E(n,:)=0;
            E(n,I)=1;
            E(I,:)=0;
            E(I,n)=1;
        end
        
        
        
        function x = Gauss_Jordan(obj, A, b)
            A = [A b];
            [m, n] = size(A);
            x = zeros(n-1,1);
            for j = 1 : m-1
                for z = 2 : m
                    if A(j,j) == 0
                        dummy = A(1,:);
                        A(1,:) = A(z,:);
                        A(z,:) = dummy;
                    end
                end
                for i = j+1:m
                    factor = A(i,j) / A(j,j);
                    A(i,:) = A(i,:) - (factor * A(j,:));
                end
            end
            for j = m : -1 : 2
                for i = j-1: -1 : 1
                    factor = A(i,j) / A(j,j);
                    A(i,:) = A(i,:) - (factor * A(j,:));
                end
            end
            for s = 1:m
                A(s,:) = A(s,:)/A(s,s);
                x(s) = A(s,n);
            end
        end
        
        function x = Gauss(obj, A, b, tol, er)
            n = length(b);
            s = zeros(n,1);
            for i = 1:n
                s(i) = abs(A(i,1));
                for j = 2:n
                    if abs(A(i,j)) > s(i)
                        s(i) = abs(A(i,j));
                    end
                end
            end
            [A, b, er] = obj.eliminate(A, b, n, s, tol, er);
            if er ~= -1
                x = obj.back_sub(A, n, b);
            end
        end
        
        function [A, b, er] = eliminate(obj, A, b, n, s, tol, er)
            for k = 1 : n-1
                [A, b, s] = obj.Pivot(A, b, s, n, k);
                if abs(A(k,k)/s(k)) < tol
                    er = -1;
                    return;
                end
                
                for i = k+1:n
                    factor = A(i,k) / A(k,k);
                    for j = k+1:n
                        A(i,j) = A(i,j) - factor*A(k,j);
                    end
                    b(i) = b(i) - factor * b(k);
                end
            end
            if abs(A(n,n)/s(n)) < tol
                er = -1 ;
            end
        end
        
        function [A, b, s] = Pivot(obj, A, b, s, n, k)
            p = k;
            big = abs(A(k,k)/s(k));
            for i = k+1:n
                dummy = abs(A(i,k)/s(i));
                if dummy > big
                    big = dummy;
                    p = i;
                end
            end
            if p ~= k
                for j = k:n
                    dummy = A(p,j);
                    A(p,j) = A(k,j);
                    A(k,j) = dummy;
                end
                dummy = b(p);
                b(p) = b(k);
                b(k) = dummy;
                
                dummy = s(p);
                s(p) = b(k);
                s(k) = dummy;
            end
        end
        
        function x = back_sub(obj, A, n, b)
            x = zeros(n,1);
            x(n) = b(n) / A(n,n);
            for i = n-1 : -1 : 1
                sum = 0;
                for j = i+1 : n
                    sum = sum + (A(i,j) * x(j));
                end
                x(i) = (b(i) - sum) / A(i,i);
            end
        end
        
        function y = forward_sub(obj, L,w)
            [nRow, nCol] = size(L);
            y = zeros(nRow, 1);
            y(1) = w(1)/L(1,1);
            w=w(:);
            for n = 2:nRow
                y(n) = (w(n) - L(n,1:n-1)*y(1:n-1))/L(n,n);
            end
        end
        
        
    end
end