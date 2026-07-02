% Copyright (c) 2026 Matteo Giordano and Kolyan Ray
%
% Code accompanying the article "Semiparametric Bernstein-von Mises theorems
% for reversible diffusions"
% by Matteo Giordano and Kolyan Ray
%
% For MATLAB R2026a
%
% Requires PDE Toolbox, Wavelet Toolbox, and Craveri package for 
% scientifically accurate coloring

%%
% Contents:
% 1. Continuously-observed bidimensional reversible diffusion with periodic
%    gradient drift vector field
%
%   1.1 Domain and ground truth
%   1.2 Continuous diffusion path
%
% 2. Posterior inference on B and Psi(B) with a truncated Besov-Laplace wavelet
%    series prior
%
%   2.1 Likelihood construction in wavelet coordinates
%   2.2 Prior draw
%   2.3 MLE and MAP
%   2.4 wpCN posterior sampler
%   2.5 Monte Carlo approximation of the plug-in posterior distributions

%%
% 1. Continuously-observed bidimensional reversible diffusion with periodic
%    gradient drift vector field

%%
% Specify the domain and create mesh for function discretisation and plots

% Display more digits
%format long

% Create PDE model
model = createpde();

% Define square domain with centre (0.5,0.5) and side length 1
SQ = [3;4;0;1;1;0;0;0;1;1];

% Create geometry
geom = decsg(SQ); 
geometryFromEdges(model,geom);

% Create mesh
Hmax = 0.025;
generateMesh(model, 'Hmax', Hmax);
mesh_nodes = model.Mesh.Nodes; 
    % 2 x mesh_size matrix whose columns contain the (x,y) coordinates of 
    % the nodes in the mesh
xn = mesh_nodes(1,:); yn = mesh_nodes(2,:);
    % x and y coordinates of nodes
mesh_nodes_num = size(mesh_nodes,2); 
    % number of nodes in the mesh
mesh_elements = model.Mesh.Elements; 
    % 6 x mesh_elements_num matrix whose columns contain the 6 node indices 
    % identifying each triangle. The first 3 elements of each column 
    % contain the indices of the 3 vertices of the triangle 
mesh_elements_num = size(mesh_elements,2); 
    % number of triangles in the mesh
[~,mesh_elements_area] = area(model.Mesh);
    % area of triangles in the mesh

% Plot mesh
%figure()
%axes('FontSize', 20, 'NextPlot','add')
%pdeplot(model)
%axis equal
%xticks([0, 0.25, 0.5, 0.75, 1]);
%yticks([0, 0.25, 0.5, 0.75, 1]);
%xlabel('x', 'FontSize', 20);
%ylabel('y', 'FontSize', 20);

% Compute triangle barycenters
tri = mesh_elements(1:3, :);
barycenters = (mesh_nodes(:, tri(1,:)) + mesh_nodes(:, tri(2,:)) + mesh_nodes(:, tri(3,:))) / 3;
xb = barycenters(1,:); yb = barycenters(2,:);
    % x and y coordinates of barycenters

% Plot the triangular mesh with the barycenters
%figure()
%axes('FontSize', 20, 'NextPlot','add')
%pdeplot(model)
%axis equal
%hold on
%plot(barycenters(1,:),barycenters(2,:),'ok','Marker','*',...
% 'Color','black','LineWidth',1)
%title('Mesh elements and their barycenters','FontSize',20)
%xticks([0, 0.25, 0.5, 0.75, 1]);
%yticks([0, 0.25, 0.5, 0.75, 1]);
%xlabel('x', 'FontSize', 20);
%ylabel('y', 'FontSize', 20);

%%
% Specify periodic potential

truth_id = 1;

switch truth_id
    case 1
        B0_raw = @(x,y) exp(-(7.5*x-5).^2 - (7.5*y-5).^2) + ...
        exp(-(7.5*x-2.5).^2 - (7.5*y-2.5).^2);

        GradB0 = @(x,y) ...
        -2*exp(-(7.5*x-5).^2 - (7.5*y-5).^2) .* ...
        [7.5*(7.5*x-5), 7.5*(7.5*y-5)] ...
        -2*exp(-(7.5*x-2.5).^2 - (7.5*y-2.5).^2) .* ...
        [7.5*(7.5*x-2.5), 7.5*(7.5*y-2.5)];

    case 2
        B0_raw = @(x,y) 2 + exp(-(7.5*x-5).^2 - (7.5*y-5).^2) - ...
        exp(-(7.5*x-2.5).^2 - (7.5*y-2.5).^2);

        GradB0 = @(x,y) ...
        -2*exp(-(7.5*x-5).^2 - (7.5*y-5).^2) .* ...
        [7.5*(7.5*x-5), 7.5*(7.5*y-5)] ...
        +2*exp(-(7.5*x-2.5).^2 - (7.5*y-2.5).^2) .* ...
        [7.5*(7.5*x-2.5), 7.5*(7.5*y-2.5)];

    case 3

        B0_raw = @(x,y) exp(-(7.5*x-5.5).^2 - (7.5*y-5.5).^2) + ...
        0.75*exp(-(5*x-1.25).^2 - (7.5*y-5.5).^2) + ...
        1.25*exp(-(7.5*x-5.5).^2 - (5*y-1.25).^2) + ...
        exp(-(7.5*x-2).^2 - (7.5*y-2).^2);

        GradB0 = @(x,y) ...
        -2*exp(-(7.5*x-5.5).^2 - (7.5*y-5.5).^2) .* ...
        [7.5*(7.5*x-5.5), 7.5*(7.5*y-5.5)] ...
        -2*0.75*exp(-(5*x-1.25).^2 - (7.5*y-5.5).^2) .* ...
        [5*(5*x-1.25), 7.5*(7.5*y-5.5)] ...
        -2*1.25*exp(-(7.5*x-5.5).^2 - (5*y-1.25).^2) .* ...
        [7.5*(7.5*x-5.5), 5*(5*y-1.25)] ...
        -2*exp(-(7.5*x-2).^2 - (7.5*y-2).^2) .* ...
        [7.5*(7.5*x-2), 7.5*(7.5*y-2)];
end

% Center B0 by subtracting its numerical integral
B0_raw_bary = B0_raw(barycenters(1,:), barycenters(2,:));
B0_int = sum(B0_raw_bary .* mesh_elements_area);
B0 = @(x,y) B0_raw(x,y) - B0_int;

% Evaluate centered B0
B0_mesh = B0(mesh_nodes(1,:), mesh_nodes(2,:));
B0_bary = B0(barycenters(1,:), barycenters(2,:));
B0_int = sum(B0_bary.*mesh_elements_area);

% L2 norm of centered B0 at the mesh nodes and barycenters
B0_norm = sqrt(sum(B0_bary.^2 .* mesh_elements_area));

% Plot B_0
clim_vals = [min(B0_mesh), max(B0_mesh)];
figure()
axes('FontSize', 20, 'NextPlot','add')
pdeplot(model,'XYData',real(B0_mesh),'ColorMap',jet)
axis equal
% pdeplot(model,'XYData',B0_mesh,'ZData',B0_mesh,'ColorMap',jet) 
    % 3D plot
title('True potential B_0', 'FontSize', 20);
xticks([0, 0.25, 0.5, 0.75, 1]);
yticks([0, 0.25, 0.5, 0.75, 1]);
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
crameri vik
cb = colorbar;
cb.FontSize = 20;
cb.Ticks = linspace(clim_vals(1), clim_vals(2), 5); % fewer ticks
cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false);

%%
% Projection of B0 onto a periodic wavelet basis

% Wavelet settings
Jgrid = 4;
N = 2^Jgrid;
    % dyadic grid size: N x N, with N = 2^Jgrid
Jwave = Jgrid;              
    % retained wavelet decomposition depth
wname = 'db6';          
    % Daubechies wavelets
dwtmode('per','nodisp');
    % periodic boundary handling for the DWT

% Create dyadic tensor grid on T^2 = [0,1)^2
xw = (0:N-1)/N;
[yw1, yw2] = meshgrid(xw, xw);
xw_ext = [xw, 1];
xn_per = mod(xn, 1);
yn_per = mod(yn, 1);
xb_per = mod(xb, 1);
yb_per = mod(yb, 1);

% Evaluate B0 on the NXN wavelet grid
B0_wave_grid = B0(yw1, yw2);
B0_wave_grid = B0_wave_grid - mean(B0_wave_grid(:));
    % re-center on the wavelet grid to remove the constant coefficient

% Periodic 2D wavelet decomposition
[C0_wave, S0_wave] = wavedec2(B0_wave_grid, Jwave, wname);
    % C0_wave is a 1 x N^2 row vector containing the wavelet coefficients
    % in the discrete periodic wavelet basis on the N x N grid

% Reconstruction from the retained wavelet coefficients
B0_wave_proj_grid = waverec2(C0_wave, S0_wave, wname);
    % NxN matrix containing the values of the wavelet projection at the NxN
    % wavelet grid
B0_wave_proj_grid = B0_wave_proj_grid - mean(B0_wave_proj_grid(:));
    % numerical re-centering

% Interpolation of wavelet projection back to mesh nodes and barycenters
B0_wave_interp = B0_wave_proj_grid';
B0_wave_ext = [B0_wave_interp, B0_wave_interp(:,1)];
B0_wave_ext = [B0_wave_ext; B0_wave_ext(1,:)];
F_wave = griddedInterpolant({xw_ext, xw_ext}, B0_wave_ext, 'spline', 'nearest');
B0_wave_mesh = F_wave(xn_per, yn_per);
B0_wave_bary = F_wave(xb_per, yb_per);
B0_wave_int = sum(B0_wave_bary .* mesh_elements_area);
B0_wave_mesh = B0_wave_mesh - B0_wave_int;
B0_wave_bary = B0_wave_bary - B0_wave_int;

% Plot B0 and projection
figure();
subplot(1,2,1)
pdeplot(model, 'XYData', real(B0_mesh));
axis equal
title('True potential B_0', 'FontSize', 20);
xlim([0 1]); ylim([0 1]);
xticks([0 0.25 0.5 0.75 1]);
yticks([0 0.25 0.5 0.75 1]);
xlabel('x', 'FontSize', 10);
ylabel('y', 'FontSize', 10);
colormap(crameri('vik'))
cb = colorbar;
cb.FontSize = 10;
cb.Ticks = linspace(clim_vals(1), clim_vals(2), 5);
cb.TickLabels = arrayfun(@(z) sprintf('%.2f', z), cb.Ticks, 'UniformOutput', false);

subplot(1,2,2)
pdeplot(model, 'XYData', B0_wave_mesh);
axis equal
title('Projection of B_0', 'FontSize', 20);
xlim([0 1]); ylim([0 1]);
xticks([0 0.25 0.5 0.75 1]);
yticks([0 0.25 0.5 0.75 1]);
xlabel('x', 'FontSize', 10);
ylabel('y', 'FontSize', 10);
colormap(crameri('vik'))
cb = colorbar;
cb.FontSize = 10;
cb.Ticks = linspace(clim_vals(1), clim_vals(2), 5);
cb.TickLabels = arrayfun(@(z) sprintf('%.2f', z), cb.Ticks, 'UniformOutput', false);

approx_error = sqrt(sum((B0_bary - B0_wave_bary).^2 .* mesh_elements_area));
disp(['L^2 approximation error via projection = ', num2str(approx_error)]);
disp(['Relative approximation error via projection = ', num2str(approx_error / B0_norm)]);

% Build index map for packed wavelet coefficient vector
wave_index = table();
start_idx = 1;

% Approximation block A_Jwave
A_size = S0_wave(1,:);
numA = prod(A_size);
    % size of coarse approximation

idx = start_idx:(start_idx + numA - 1);
wave_index = [wave_index;

table("A", Jwave, A_size(1), A_size(2), ...
    idx(1), idx(end), ...
    'VariableNames', {'Type','Level','Rows','Cols','Start','End'})];

start_idx = start_idx + numA;

% Detail blocks: H, V, D
for row = 2:size(S0_wave,1)-1
    level = Jwave - row + 2;
    block_size = S0_wave(row,:);
    numBlock = prod(block_size);
    for typ = ["H", "V", "D"]
        idx = start_idx:(start_idx + numBlock - 1);
        wave_index = [wave_index;
        table(typ, level, block_size(1), block_size(2), ...
        idx(1), idx(end), ...
        'VariableNames', {'Type','Level','Rows','Cols','Start','End'})];
        start_idx = start_idx + numBlock;
    end
end

% Display wavelet index structure
% disp(wave_index)

% Identify active coefficients (remove constant mode)
num_wave_coeffs = length(C0_wave);
constant_idx = wave_index.Start( ...
    wave_index.Type == "A" & wave_index.Level == Jwave);

active_idx = setdiff((1:num_wave_coeffs)', constant_idx);
num_active_coeffs = length(active_idx);

%%
% 1.2 Continuous diffusion path

%%
% Sample path simulation

rng(100);

% Generate Euler–Maruyama approximations
sigma = 1; 
% standard deviation of Brownian motion
sim_len = 200000; 
% number of simulated values
deltaT = .0005;
% small time interval between Euler-Maruyama iterates
T = sim_len *deltaT;
% time horizon

% Initialise path X_0, ..., X_sim_len
X = zeros(2,sim_len + 1); 
x_0 = [.5,.5]; X(:,1) = x_0;

% Pre-generate Brownian increments
dW = sqrt(deltaT) * randn(2, sim_len);

for i=1:sim_len

    % Project current point onto torus for periodisation
    X_proj = mod(X(:,i), 1);

    % Euler-Maruyama step
    grad = GradB0(X_proj(1), X_proj(2))';
    X(:,i+1) = X(:,i) + grad * deltaT + sigma * dW(:,i);
end

% Plot trajectory
figure()
axes('FontSize', 20, 'NextPlot','add')
plot(X(1,:),X(2,:),'blue','LineWidth',1);
plot(X(1,1),X(2,1),'Marker','*','Color','red','LineWidth',2)
plot(X(1,end),X(2,end),'Marker','*','Color','yellow','LineWidth',2)
title(['T=',num2str(T)],'Fontsize',20)
legend('X_t','X_0','X_T','Fontsize',20)
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);

% Plot trajectory with gradient vector field
%Compute gradient
%[xg, yg] = meshgrid(linspace(min(X(1,:)),max(X(1,:)),100), linspace(min(X(2,:)),max(X(2,:)),100));
%Gx = zeros(size(xg));
%Gy = zeros(size(yg));
%for j = 1:numel(xg)
%    g = GradB0(mod(xg(j),1), mod(yg(j),1));
%    Gx(j) = g(1);
%    Gy(j) = g(2);
%end
%figure();
%axes('FontSize', 20, 'NextPlot', 'add');
%plot(X(1,:), X(2,:), 'blue', 'LineWidth', 1);
%plot(X(1,1), X(2,1),'Marker', '*', 'Color', 'red', 'LineWidth', 2);
%plot(X(1,end), X(2,end),'Marker', '*', 'Color', 'yellow', 'LineWidth', 2);
%quiver(xg, yg, Gx, Gy, 'LineWidth', 1, 'Color', 'green');
%axis equal;
%title(['T=',num2str(T)],'Fontsize',20)
%legend('X_t', 'X_0', 'X_T', '\nabla B_0', 'FontSize', 20);
%xlabel('x', 'FontSize', 20);
%ylabel('y', 'FontSize', 20);

%%
% Compute matrix Sigma and vector H for likelihood evaluation

% Initialise gradient values along diffusion path
Gx = zeros(sim_len, num_active_coeffs);
Gy = zeros(sim_len, num_active_coeffs);
dx = 1 / N;
X1_path = mod(X(1,1:sim_len)', 1);
X2_path = mod(X(2,1:sim_len)', 1);

for k = 1:num_active_coeffs
    global_idx = active_idx(k);

    e_wave = zeros(size(C0_wave));
    e_wave(global_idx) = 1;

    phi_grid = waverec2(e_wave, S0_wave, wname);

    dphi_dx = (circshift(phi_grid, [0, -1]) - circshift(phi_grid, [0, 1])) / (2 * dx);
    dphi_dy = (circshift(phi_grid, [-1, 0]) - circshift(phi_grid, [1, 0])) / (2 * dx);

    dphi_dx_interp = dphi_dx';
    dphi_dy_interp = dphi_dy';

    dphi_dx_ext = [dphi_dx_interp, dphi_dx_interp(:,1)];
    dphi_dx_ext = [dphi_dx_ext; dphi_dx_ext(1,:)];

    dphi_dy_ext = [dphi_dy_interp, dphi_dy_interp(:,1)];
    dphi_dy_ext = [dphi_dy_ext; dphi_dy_ext(1,:)];

    Fx = griddedInterpolant({xw_ext, xw_ext}, dphi_dx_ext, 'linear', 'nearest');
    Fy = griddedInterpolant({xw_ext, xw_ext}, dphi_dy_ext, 'linear', 'nearest');

    Gx(:,k) = Fx(X1_path, X2_path);
    Gy(:,k) = Fy(X1_path, X2_path);
end

% Approximate Sigma
Sigma_wave = deltaT * (Gx' * Gx + Gy' * Gy);
Sigma_active = 0.5 * (Sigma_wave + Sigma_wave');

% Trajectory increments
dX1 = diff(X(1,1:sim_len+1))';
dX2 = diff(X(2,1:sim_len+1))';

% Stochastic integral approximation
H_wave = Gx' * dX1 + Gy' * dX2;
H_active = H_wave(:);

% wavelet coefficients-to-loglikelihood function
loglik_coeffs = @(C) H_active' * C(:) - 0.5 * (C(:)' * Sigma_active * C(:));

C0_active = C0_wave(active_idx);
C0_active = C0_active(:);
loglik_B0 = loglik_coeffs(C0_active);
disp(['Log-likelihood of projection of B_0 = ', num2str(loglik_B0)]);

%%
% 2. Posterior inference on B and Psi(B) with a truncated Besov-Laplace wavelet
%    series prior

%%
% 2.1 Besov-Laplace prior

d = 2;
    % dimension
s = 1.25;
    % smoothness parameter
%prior_global_scale = T^(-d/(2*s + d));
prior_global_scale = 1^(-d/(2*s + d));

% Build scale vector tau_j = T^{-d/(2s+d)} 2^{-l(j)(s+1-d/2)} = T^{-2/(2s+2)} 
% 2^{-l(j)s}
tau_wave = zeros(size(C0_wave));

for b = 1:height(wave_index)
    idx = wave_index.Start(b):wave_index.End(b);

    if wave_index.Type(b) == "A"
        lev = 0;
    else
        % MATLAB detail level 1 is finest; paper level l large is finest.
        lev = Jwave - wave_index.Level(b) + 1;
    end

    tau_wave(idx) = prior_global_scale * 2^(-lev * (s + 1 - d/2));
end


% Remove coefficients for constants
tau_wave(constant_idx) = 0;
tau_active = tau_wave(active_idx);
tau_active = tau_active(:);
lambda_active = 1 ./ tau_active;

% Draw one Laplace prior sample
U = min(max(rand(size(C0_wave)), eps), 1-eps) - 0.5;
G = -sign(U) .* log(1 - 2*abs(U));
    % transformation to Laplace random variables
C_prior = tau_wave .* G;
C_prior(constant_idx) = 0;

% Reconstruct prior draw on mesh and barycenters via coeffs_to_mesh() function
[B_prior_mesh, B_prior_bary, ~] = coeffs_to_mesh(C_prior(active_idx), active_idx, ...
    C0_wave, S0_wave, wname, xw_ext, xn_per, yn_per, xb_per, yb_per, mesh_elements_area);

% log-likelihood of prior draw
loglik_prior_draw = loglik_coeffs(C_prior(active_idx));
disp(['Log-likelihood of prior draw = ', num2str(loglik_prior_draw)]);

% Plot prior draw
figure()
axes('FontSize', 20, 'NextPlot', 'add')
pdeplot(model, 'XYData', real(B_prior_mesh), 'ColorMap', jet)
axis equal
title('Besov-Laplace prior draw', 'FontSize', 20)
xticks([0, 0.25, 0.5, 0.75, 1]);
yticks([0, 0.25, 0.5, 0.75, 1]);
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
colormap(crameri('vik'))
cb = colorbar;
clim_prior = [min(B_prior_mesh), max(B_prior_mesh)];
cb.FontSize = 20;
cb.Ticks = linspace(clim_prior(1), clim_prior(2), 5);
cb.TickLabels = arrayfun(@(z) sprintf('%.2f', z), cb.Ticks, 'UniformOutput', false);

%%
% 2.2 MLE and MAP

% MLE via explicit solution
ridge = 1e-10;
C_mle = (Sigma_active + ridge * eye(size(Sigma_active))) \ H_active;
loglik_mle = loglik_coeffs(C_mle);
disp(['log-likelihood of MLE = ', num2str(loglik_mle)]);

% Reconstruct MLE on mesh and barycenters via coeffs_to_mesh() function
[B_mle_mesh, B_mle_bary, ~] = coeffs_to_mesh(C_mle, active_idx, C0_wave, ...
    S0_wave, wname, xw_ext, xn_per, yn_per, xb_per, yb_per, mesh_elements_area);

mle_error = sqrt(sum((B_mle_bary - B0_bary).^2 .* mesh_elements_area));
disp(['MLE L^2 estimation error = ', num2str(mle_error)]);
disp(['MLE relative error = ', num2str(mle_error / B0_norm)]);

figure()
subplot(1,2,1)
pdeplot(model, 'XYData', real(B0_mesh))
axis equal
title('True potential B_0')
clim(clim_vals)
xlim([0 1]); ylim([0 1]);
xticks([0 0.25 0.5 0.75 1]);
yticks([0 0.25 0.5 0.75 1]);
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
colormap(crameri('vik'))
cb = colorbar;
cb.FontSize = 20;

subplot(1,2,2)
pdeplot(model, 'XYData', real(B_mle_mesh))
axis equal
title('MLE')
clim(clim_vals)
xlim([0 1]); ylim([0 1]);
xticks([0 0.25 0.5 0.75 1]);
yticks([0 0.25 0.5 0.75 1]);
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
colormap(crameri('vik'))
cb = colorbar;
cb.FontSize = 20;

% MAP via FISTA
max_iter_map = 10000;
tol_map = 1e-8;
    % FISTA parameters
soft_thresh = @(x,a) sign(x) .* max(abs(x) - a, 0);
    % soft-thresholding operator
L_fista = eigs(Sigma_active, 1, 'largestreal');
    % Lipschitz constant of gradient of smooth part

% FISTA initialisation. Use cold start (zero or prior draw) or MLE. Use C_0
% for debugging
C_map = zeros(num_active_coeffs,1);
    % cold start
% C_map = C0_active;
    % warm start
% C_map = C_mle;

Y = C_map;
t_fista = 1;
obj_trace_map = zeros(max_iter_map,1);

for k = 1:max_iter_map
    C_old = C_map;

    % Gradient of smooth negative log-likelihood
    grad_Y = Sigma_active * Y - H_active;

    % Proximal weighted L1 update
    C_map = soft_thresh(Y - grad_Y / L_fista, lambda_active / L_fista);

    % FISTA acceleration
    t_new = (1 + sqrt(1 + 4*t_fista^2)) / 2;
    Y = C_map + ((t_fista - 1) / t_new) * (C_map - C_old);
    t_fista = t_new;

    % Objective value
    obj_trace_map(k) = 0.5 * C_map' * Sigma_active * C_map ...
        - H_active' * C_map + sum(lambda_active .* abs(C_map));

    % Stopping rule
    rel_change = norm(C_map - C_old) / max(1, norm(C_old));
    if rel_change < tol_map
        obj_trace_map = obj_trace_map(1:k);
        disp(['FISTA converged after ', num2str(k), ' iterations.'])
        break
    end

    if k == max_iter_map
        disp('FISTA reached maximum number of iterations.')
    end
end

% FISTA diagnostics
%figure()
%plot(obj_trace_map, 'LineWidth', 1)
%xlabel('Iteration')
%ylabel('Objective')
%title('FISTA objective trace for Besov-Laplace MAP')
%grid on

loglik_map = loglik_coeffs(C_map);
logprior_map = -sum(lambda_active .* abs(C_map));
logpost_map = loglik_map + logprior_map;
disp(['Log-likelihood of MAP = ', num2str(loglik_map)])
%disp(['MAP log-prior penalty = ', num2str(logprior_map)])
%disp(['MAP log-posterior = ', num2str(logpost_map)])

% Reconstruct MAP on mesh and barycenters via coeffs_to_mesh() function
[B_map_mesh, B_map_bary, ~] = coeffs_to_mesh(C_map, active_idx, C0_wave, ...
    S0_wave, wname, xw_ext, xn_per, yn_per, xb_per, yb_per, mesh_elements_area);

map_error = sqrt(sum((B_map_bary - B0_bary).^2 .* mesh_elements_area));
disp(['MAP L^2 estimation error = ', num2str(map_error)]);
disp(['MAP relative error = ', num2str(map_error / B0_norm)]);

% Plot MAP
figure()
subplot(1,2,1)
pdeplot(model, 'XYData', real(B0_mesh))
axis equal
title('True potential B_0')
clim(clim_vals)
xlim([0 1]); ylim([0 1]);
xticks([0 0.25 0.5 0.75 1]);
yticks([0 0.25 0.5 0.75 1]);
xlabel('x', 'FontSize', 10);
ylabel('y', 'FontSize', 10);
colormap(crameri('vik'))
cb = colorbar;
cb.FontSize = 10;

subplot(1,2,2)
pdeplot(model, 'XYData', real(B_map_mesh))
axis equal
title('Besov-Laplace MAP estimator')
clim(clim_vals)
xlim([0 1]); ylim([0 1]);
xticks([0 0.25 0.5 0.75 1]);
yticks([0 0.25 0.5 0.75 1]);
xlabel('x', 'FontSize', 10);
ylabel('y', 'FontSize', 10);
colormap(crameri('vik'))
cb = colorbar;
cb.FontSize = 10;

%%
% 2.3 Posterior inference via whitened pCN

% Whitening transformation for wpCN
white_to_coeffs = @(w) tau_active .* local_white_to_laplace(w);

% MCMC parameters
mcmc_steps = 100000;
burnin = 10000;
thin = 10;
num_saved = floor((mcmc_steps - burnin)/thin);

% Adaptive wpCN parameters
target_accept = 0.30;
b_wpCN = 0.05;
        % initial stepsize
b_min = 1e-5; b_max = 0.49;
        % range for stepsize
adapt_interval = 50;
adapt_strength = 0.5;

% Initialisation
w_current = zeros(num_active_coeffs,1);
    % cold start

% Warm start at MAP
%g0 = C_map ./ tau_active;
%p0 = zeros(size(g0));
%pos = g0 >= 0;
%p0(pos) = 1 - 0.5*exp(-g0(pos));
%p0(~pos) = 0.5*exp(g0(~pos));
%p0 = min(max(p0, eps), 1-eps);
%w_current = norminv(p0);

C_current = white_to_coeffs(w_current);
loglik_current = loglik_coeffs(C_current);

% Storage for posterior samples
C_samples = zeros(num_active_coeffs, num_saved);

% Storage for diagnostic
loglik_samples = zeros(num_saved,1);
loglik_trace = zeros(mcmc_steps,1);
accept_trace = zeros(mcmc_steps,1);
accept_step = zeros(mcmc_steps,1);
b_trace = zeros(mcmc_steps,1);

accept_count = 0;
accept_block = 0;
save_counter = 0;

rng(500);

for it = 1:mcmc_steps

    % Whitened pCN proposal
    zeta = randn(num_active_coeffs,1);
    w_prop = sqrt(1 - 2*b_wpCN) * w_current + sqrt(2*b_wpCN) * zeta;

    % Transform proposal to Besov-Laplace wavelet coefficients
    C_prop = white_to_coeffs(w_prop);
    loglik_prop = loglik_coeffs(C_prop);

    % wpCN acceptance probability via likelihood ratio
    log_alpha = loglik_prop - loglik_current;
    accepted = 0;

    if log(rand) < log_alpha
        w_current = w_prop;
        C_current = C_prop;
        loglik_current = loglik_prop;
        accept_count = accept_count + 1;
        accepted = 1;
    end

    % Adapt step size during burn-in only
    accept_block = accept_block + accepted;

    if it <= burnin && mod(it, adapt_interval) == 0
        block_accept = accept_block / adapt_interval;

        % Update b_wpCN on logit scale, constrained to (0,b_max)
        logit_b = log(b_wpCN / (b_max - b_wpCN));
        logit_b = logit_b + adapt_strength * (block_accept - target_accept);
        b_wpCN = b_max * exp(logit_b) / (1 + exp(logit_b));
        b_wpCN = min(max(b_wpCN, b_min), b_max);
        accept_block = 0;
    end

    % Register diagnostics
    loglik_trace(it) = loglik_current;
    accept_trace(it) = accept_count / it;
    accept_step(it) = accepted;
    b_trace(it) = b_wpCN;

    % Save thinned posterior samples after burn-in
    if it > burnin && mod(it - burnin, thin) == 0
        save_counter = save_counter + 1;
        C_samples(:,save_counter) = C_current;
        loglik_samples(save_counter) = loglik_current;
    end
end

% Acceptance rate
accept_rate = accept_count / mcmc_steps;
post_burn_accept_rate = mean(accept_step(burnin+1:end));
%disp(['Overall wpCN acceptance rate = ', num2str(accept_rate)]);
disp(['Post-burn-in wpCN acceptance rate = ', num2str(post_burn_accept_rate)]);
%disp(['Final b_wpCN = ', num2str(b_wpCN)]);

%%
% Plot of wpCN diagnostics and posterior mean

% wpCN diagnostics
figure()
subplot(3,1,1)
plot(1:mcmc_steps, loglik_trace, 'LineWidth', 1)
hold on
yline(loglik_B0, 'r', 'LineWidth', 2)
xlabel('Iteration')
ylabel('Log-likelihood')
title('wpCN log-likelihood trace')
legend({'$\ell_T(B_n)$', '$\ell_T(B_0)$'}, ...
    'Interpreter', 'latex', ...
    'Location', 'best')
grid on

subplot(3,1,2)
plot(1:mcmc_steps, accept_trace, 'LineWidth', 1)
hold on
yline(target_accept, 'g', 'LineWidth', 2)
xlabel('Iteration')
ylabel('Cumulative acceptance rate')
title('wpCN cumulative acceptance rate')
ylim([0 1])
legend('Acceptance rate', 'Target 0.30', 'Location', 'best')
grid on

subplot(3,1,3)
plot(1:mcmc_steps, b_trace, 'LineWidth', 1)
xlabel('Iteration')
ylabel('b_{wpCN}')
title('Adaptive wpCN step size')
grid on

% Posterior mean
post_mean_wave = mean(C_samples, 2);
[post_mean_mesh, post_mean_bary, ~] = coeffs_to_mesh(post_mean_wave, active_idx, ...
    C0_wave, S0_wave, wname, xw_ext, xn_per, yn_per, xb_per, yb_per, mesh_elements_area);

% Plot posterior mean
figure()
subplot(1,2,1)
pdeplot(model, 'XYData', real(B0_mesh))
axis equal
title('True potential B_0', 'FontSize', 20);
clim(clim_vals)
xlim([0 1]); ylim([0 1]);
xticks([0 0.25 0.5 0.75 1]);
yticks([0 0.25 0.5 0.75 1]);
xlabel('x', 'FontSize', 10);
ylabel('y', 'FontSize', 10);
colormap(crameri('vik'))
cb = colorbar;
cb.FontSize = 10;

subplot(1,2,2)
pdeplot(model, 'XYData', post_mean_mesh)
axis equal
title('Posterior mean', 'FontSize', 20);
xlim([0 1]); ylim([0 1]);
xticks([0 0.25 0.5 0.75 1]);
yticks([0 0.25 0.5 0.75 1]);
xlabel('x', 'FontSize', 10);
ylabel('y', 'FontSize', 10);
colormap(crameri('vik'));
cb = colorbar;
cb.FontSize = 10;
cb.Ticks = linspace(clim_vals(1), clim_vals(2), 5);
cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false);
clim(clim_vals);

est_error = sqrt(sum((B0_bary - post_mean_bary).^2 .* mesh_elements_area));
disp(['Posterior mean L^2 estimation error = ', num2str(est_error)]);
disp(['Posterior mean relative estimation error = ', num2str(est_error / B0_norm)]);
disp(['L^2 approximation error via projection = ', num2str(approx_error)]);
disp(['Relative approximation error via projection = ', num2str(approx_error / B0_norm)]);

%%
% 2.4 Plug-in posterior distributions for functionals

num_post_samples = size(C_samples, 2);
w = mesh_elements_area(:);

% True functionals
Psi1_true = sum((B0_bary(:).^2) .* w);
q = 4;
Psi2_true = sum((B0_bary(:).^q) .* w);

exp2B0 = exp(1 * B0_bary(:));
Z_B0 = sum(exp2B0 .* w);
mu_B0 = exp2B0 / Z_B0;
Psi3_true = sum(mu_B0 .* log(mu_B0) .* w);

Psi_true = [Psi1_true; Psi2_true; Psi3_true];

% Precompute wavelet basis matrix at barycenters
Phi_wave_bary = zeros(mesh_elements_num, num_active_coeffs);

xlims = [0.02 0.07; 0.00 0.03; 0.01 0.045];

for k = 1:num_active_coeffs
    e_wave = zeros(size(C0_wave));
    e_wave(active_idx(k)) = 1;

    phi_grid = waverec2(e_wave, S0_wave, wname);

    phi_interp = phi_grid';
    phi_ext = [phi_interp, phi_interp(:,1)];
    phi_ext = [phi_ext; phi_ext(1,:)];

    F_phi = griddedInterpolant({xw_ext, xw_ext}, phi_ext, ...
        'spline', 'nearest');

    Phi_wave_bary(:,k) = F_phi(xb_per, yb_per);
end

% Evaluate all posterior draws at barycenters
B_draws_bary = Phi_wave_bary * C_samples;

% Plug-in posterior samples
Psi_samples = zeros(num_post_samples, 3);

Psi_samples(:,1) = sum((B_draws_bary.^2) .* w, 1)';
Psi_samples(:,2) = sum((B_draws_bary.^q) .* w, 1)';

exp2B = exp(1 * B_draws_bary);
Z_B = sum(exp2B .* w, 1);
mu_B = exp2B ./ Z_B;

Psi_samples(:,3) = sum(mu_B .* log(mu_B) .* w, 1)';

% Posterior summaries
Psi_post_mean = mean(Psi_samples, 1);
Psi_post_sd = std(Psi_samples, 0, 1);
Psi_CI = quantile(Psi_samples, [0.025 0.975], 1);

functional_titles = { ...
    '\Psi(B) = \int_{[0,1]^d} B^2(x)dx', ...
    '\Psi(B) = \int_{[0,1]^d} B^4(x)dx', ...
    '\Psi(B) = \int_{[0,1]^d} log[\mu_B(x)] \mu_B(x) dx'};

for j = 1:3
    
    figure()
    axes('FontSize', 20, 'NextPlot','add')
    histogram(Psi_samples(:,j), ...
    'Normalization', 'pdf', ...
    'BinLimits', xlims(j,:), ...
    'NumBins', 40)
    xline(Psi_true(j), 'r', 'LineWidth', 2)
    xline(Psi_CI(1,j), 'b', 'LineWidth', 2)
    xline(Psi_CI(2,j), 'b', 'LineWidth', 2)
    xx = linspace(xlims(j,1), xlims(j,2), 1000);
    plot(xx, ...
    normpdf(xx, Psi_post_mean(j), Psi_post_sd(j)), ...
    'green', 'LineWidth', 2)
    xlim(xlims(j,:))
    legend({'$\Pi(\Psi(B)\mid X^T)$', '$\Psi(B_0)$'}, ...
    'Interpreter', 'latex', ...
    'FontSize', 20, ...
    'Location', 'best')
    grid on
end

%%
% Local helper functions

% Robust whitening transformation
function g = local_white_to_laplace(w)
    w = w(:);
    p = normcdf(abs(w));
    p = min(max(p, eps), 1 - eps);
    g = sign(w) .* (-log(2 - 2*p));
    g(w == 0) = 0;
end

% Reconstruction from wavelet coefficients to PDE mesh and barycenters
function [B_mesh, B_bary, B_grid] = coeffs_to_mesh(C_active, active_idx, C_template, ...
    S_wave, wname, xw_ext, xn_per, yn_per, xb_per, yb_per, mesh_elements_area)

    C_full = zeros(size(C_template));
    C_full(active_idx) = C_active(:);

    B_grid = waverec2(C_full, S_wave, wname);
    B_grid = B_grid - mean(B_grid(:));

    B_interp = B_grid';
    B_ext = [B_interp, B_interp(:,1)];
    B_ext = [B_ext; B_ext(1,:)];

    F = griddedInterpolant({xw_ext, xw_ext}, B_ext, 'spline', 'nearest');

    B_mesh = F(xn_per, yn_per);
    B_bary = F(xb_per, yb_per);

    B_int = sum(B_bary .* mesh_elements_area);
    B_mesh = B_mesh - B_int;
    B_bary = B_bary - B_int;
end
