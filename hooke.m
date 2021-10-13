clear all;
clc;
% simulation's parameters
n = 4;
k = 1;
dist = 4;
% temporal horizon
H = 100;
% step size
step = 0.1;
T = [0:step:H];
equilibrium = [0:n-1]'*dist;
% initial conditions
q0 = [1;zeros(n-2,1);-1];
x0 = equilibrium + q0;
v0 = zeros(n,1);
% exact solution
x = hooke_exact_solution(n, k, dist, H, step, T, x0, v0);
% approx solution
x_euler = hooke_euler_solution(n, k, dist, H, step, T, x0, v0);
x_verlet = hooke_verlet_solution(n, k, dist, H, step, T, x0, v0);

%particles' film
figure;
plot([],[]);
for t = [1:length(T)]
    scatter(x(:,t), zeros(n,1), 'filled');    
    xlim([-dist dist*n+dist]);
    pause(0.1);
    drawnow;
end

% plotting motion of particles
figure;
subplot(1,3,1);
plot(T, x);
title("Hooke, exact solution");
subtitle("1D, "+n+" particles");
xlabel("time");
ylabel("space");

subplot(1,3,2);
plot(T, x_euler);
title("Hooke, approx solution");
subtitle("1D, "+n+" particles");
xlabel("time");
ylabel("space");

subplot(1,3,3);
plot(T, x_verlet);
title("Hooke, approx solution");
subtitle("1D, "+n+" particles");
xlabel("time");
ylabel("space");

savefig("Figures/comparison.fig");
% plotting errors
figure;
subplot(1,2,1);
plot(T, abs(mean(x_euler) - mean(x)));
title("Hooke, Euler mean error");
subtitle("1D, "+n+" particles");
xlabel("time");
ylabel("error");

subplot(1,2,2);
plot(T, abs(mean(x_verlet) - mean(x)));
title("Hooke, Verlet mean error");
subtitle("1D, "+n+" particles");
xlabel("time");
ylabel("error");


