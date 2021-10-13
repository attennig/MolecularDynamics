function x = hooke_euler_solution(n, k, dist, H, step, T, x0, v0)
    % coefficients matrix
    M = zeros(n,n) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1) + 2*diag(ones(n,1));
    M(1,1) = 1;
    M(n,n) = 1;
    equilibrium = [0:n-1]'*dist;
    x = zeros(n,length(T));
    x(:,1) = x0;
    q = zeros(n,length(T));
    q(:,1) = x0 - equilibrium;
    v = zeros(n,length(T));
    v(:,1) = v0;
    NSTEP = H/step;
    for t = [1:NSTEP]
        x(:,t+1) = x(:,t) + step * v(:,t);
        force = - k * (M * q(:,t));
        v(:,t+1) = v(:,t) + step * force;
        q(:,t+1) = x(:,t+1) - equilibrium;    
    end
 end