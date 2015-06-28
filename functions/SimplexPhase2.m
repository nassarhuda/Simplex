%% phase 2 simplex
function [iters, z, xopt, globalSingular] = ...
    SimplexPhase2(A,b,chat,N_ids,B_ids)

%pre set certain values
globalSingular = 0;
isVZero = false;
iters = 0;
xopt = zeros(size(chat,1),1);
z = 0;
maxit = 1000;
N = A(:,N_ids);
B = A(:,B_ids);
while(~isVZero && iters<=maxit)
    
    [L U P] = lu(B);
    c_N = chat(N_ids);
    c_B = chat(B_ids);
    % Solve Bx_B = b:
    xopt(B_ids,:) = U\(L\(P*b));
    x_B = xopt(B_ids,:);
    % Solve B^Tl = c_B
    l_vec = P'*(L'\(U'\c_B));
    
    w2hat = 0;
    getBversion2;
    if ~isVZero
        %Update the indices
        N_ids(id_val) = p;
        B_ids(l_val) = q;
        
        t = xopt(p);
        xopt(p) = xopt(q);
        xopt(q) = t;
    end
    z = c_B'*x_B;
    N = A(:,N_ids);
    B = A(:,B_ids);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iters = iters + 1;
    fprintf('iteration number %d --- z value = %f\n',iters,z);
end

    function getBversion2
        % this one based on LU decomposition
        
        isSingular = true;
        sindex = 1;
        s_N = c_N - N'*l_vec;
        [s_N_sorted s_N_sortedIDS] = sort(s_N);
        
        while isSingular
            B_ids_here = B_ids;
            v = s_N_sorted(sindex);
            if v<-.00001
                id_val = s_N_sortedIDS(sindex);
                % entering index
                q = N_ids(id_val);
                % solve Bd = a_q
                a_q = A(:,q);
                d = U\(L\(P*a_q));
                m = inf;
                l_val = 0;
                for i = 1:numel(d)
                    if d(i) > 0
                        t = x_B(i)/d(i);
                        if t<m
                            m = t;
                            l_val = i;
                        end
                    end
                end
                
                % leaving index is:
                p = B_ids_here(l_val);
                LUupdatev2;
                if abs(w2hat) < 1e-9
                    sindex = sindex + 1;
                    globalSingular = globalSingular + 1;
                else
                    isSingular = false;
                end
                
            else
                isSingular = false;
                isVZero = true;
            end
        end
        
    end
    function LUupdatev2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        a_p = A(:,p);
        a_q = A(:,q);
        e_p = zeros(size(B,1),1);
        e_p(l_val) = 1;
        t_p1 = a_q - a_p;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Upadate L and U:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        m = size(U,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create matrix PI (in notes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ROWS = [1:l_val-1,l_val+1:m,l_val]';
        
        
        Uplus = U + L\((P*t_p1)*(e_p'));
        
        % instead of multiplying by permutation matrix
        Utemp = Uplus(:,ROWS);
        UU = Utemp(ROWS,:);
        
        L_lastrow = zeros(1,m);
        
        for idl = l_val:m-1
            L_lastrow(idl) = (UU(end,idl)-((L_lastrow(1:idl-1)*...
                UU(1:idl-1,idl))))/UU(idl,idl);
        end
        
        if l_val < m
            wtemp = dot(L_lastrow(l_val:m-1),UU(l_val:m-1,end));
        else wtemp = 0;
        end
        w2hat = UU(end,end) - wtemp;
    end


end



