%% Input values:
function [iters z xopt B_ids] = ...
    SimplexPhase1(A,b,chat,B_ids,N_ids,limit)

%pre set certain values
isVZero = false;
iters = 0;
xopt = zeros(size(chat,1),1);
z = 0;
maxit = 1000;
N = A(:,N_ids);
B = A(:,B_ids);

while(~isVZero && iters<=maxit)
    
    [L, U, P] = lu(B);
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
        B_ids(lval) = q;
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
                lval = 0;
                for i = 1:numel(d)
                    if d(i) > 0
                        t = x_B(i)/d(i);
                        if t<=m
                            m = t;
                            lval = i;
                        end
                    end
                end
                
                % leaving index is:
                p = B_ids_here(lval);
                LUupdatev2;
                if q > limit || abs(w2hat) < 1e-10
                    sindex = sindex + 1;
                    disp('got singular');
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
        e_p(lval) = 1;
        t_p1 = a_q - a_p;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Upadate L and U:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        m = size(U,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create matrix PI (in notes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ROWS = [1:lval-1,lval+1:m,lval]';
        
        
        Uplus = U + L\((P*t_p1)*(e_p'));
        
        %instead of multiplying by a permutation matrix
        Utemp = Uplus(:,ROWS);
        UU = Utemp(ROWS,:);
        
        L_lastrow = zeros(1,m);
        
        for idl = lval:m-1
            L_lastrow(idl) = (UU(end,idl)-((L_lastrow(1:idl-1)*...
                UU(1:idl-1,idl))))/UU(idl,idl);
        end
        
        if lval < m
            wtemp = dot(L_lastrow(lval:m-1),UU(lval:m-1,end));
        else wtemp = 0;
        end
        w2hat = UU(end,end) - wtemp;
    end

end

