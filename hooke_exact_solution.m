function x = hooke_exact_solution(n, k, dist, H, step, T, x0, v0)
    
    q0 = x0 - [0:n-1]'*dist;
    q0_dx = v0;
    % eigenvalues in matrix except the first
    E_VAL = 2*(1-cos([1:n-1]*pi/n)).*ones(n,n-1); 
    % eigenvectors except the first
    E_VEC = cos([1:n-1].*pi/n.*([1:n]'-0.5));   

    A = [E_VEC zeros(n,n-1); zeros(n,n-1) E_VEC.*sqrt(k*E_VAL)];

    % finding constants
    d = A \ [q0; q0_dx];
    c = d(1:n-1);
    C = ones(n,1)*c';
    deltas = d(n:2*(n-1));
    DELTAS = ones(n,1)*deltas';
    % solving motion equations
    q = zeros(n,length(T));
    x = zeros(n,length(T));
    for t = [1:length(T)]
        q(:,t) = sum(C.*E_VEC.*cos(t*sqrt(k*E_VAL)+DELTAS),2);
        x(:,t) = dist * [0:n-1]' + q(:,t);
    end
end