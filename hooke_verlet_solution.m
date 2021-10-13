function x = hooke_verlet_solution(n, k, dist, H, step, T, x0, v0)
    % coefficients matrix
    M = zeros(n,n) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1) + 2*diag(ones(n,1));
    M(1,1) = 1;
    M(n,n) = 1;
   
    equilibrium = [0:n-1]'*dist;
    x = zeros(n,length(T));
    x(:,1) = x0;
    q = zeros(n,length(T));
    q(:,1) = x(:,1) - equilibrium;
    
    x(:,2) = x(:,1) - step * v0 + step^2/2 * (- k * M * q(:,1));
    q(:,2) = x(:,2) - equilibrium;
    
    NSTEP = H/step;
    for t = [2:NSTEP]
        x(:,t+1) = 2*x(:,t) - x(:,t-1) + step^2*(- k * M * q(:,t));
        q(:,t+1) = x(:,t+1) - equilibrium;
    end

end