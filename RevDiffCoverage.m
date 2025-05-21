%%
% Specify the domain and create the mesh

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
generateMesh(model,'Hmax',0.025);
mesh_nodes = model.Mesh.Nodes; 
    % 2 x mesh_size matrix whose columns contain the (x,y) coordinates of 
    % the nodes in the mesh
mesh_nodes_num = size(mesh_nodes); mesh_nodes_num=mesh_nodes_num(2); 
    % number of nodes in the mesh
mesh_elements = model.Mesh.Elements; 
    % 6 x mesh_elements_num whose columns contain the 6 node indices 
    % identifying each triangle. The first 3 elements of each column 
    % contain the indices of the 3 vertices of the triangle 
mesh_elements_num = size(mesh_elements); 
mesh_elements_num = mesh_elements_num(2); 
    % number of triangles in the mesh
[~,mesh_elements_area] = area(model.Mesh);

% Plot mesh
figure()
axes('FontSize', 20, 'NextPlot','add')
pdeplot(model)
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);

% Compute barycenters of triangles
barycenters = zeros(2,mesh_elements_num);
for i=1:mesh_elements_num
    barycenters(:,i) = mean(mesh_nodes(:,mesh_elements(1:3,i)),2);
end

% Plot the triangular mesh with the barycenters
%figure()
%axes('FontSize', 15, 'NextPlot','add')
%pdeplot(model)
%hold on
%plot(barycenters(1,:),barycenters(2,:),'ok','Marker','*',...
% 'Color','black','LineWidth',1) 
%title('Mesh elements and their barycenters','FontSize',20)
%xlabel('x', 'FontSize', 15);
%ylabel('y', 'FontSize', 15);

%%
% Specify periodic potential

% Specify B_0 and its gradient as functions of (x,y)

%Ground truth 1
%B0 = @(x,y) exp(-(7.5*x-5).^2-(7.5*y-5).^2)+exp(-(7.5*x-2.5).^2-(7.5*y-2.5).^2);
%GradB0 = @(x,y) -2*exp(-(7.5*x-5).^2-(7.5*y-5).^2)*[7.5*(7.5*x-5),7.5*(7.5*y-5)]...
%    -2*exp(-(7.5*x-2.5).^2-(7.5*y-2.5).^2)*[7.5*(7.5*x-2.5),7.5*(7.5*y-2.5)];

%Ground truth 2
%B0 = @(x,y) 2+exp(-(7.5*x-5).^2-(7.5*y-5).^2)-exp(-(7.5*x-2.5).^2-(7.5*y-2.5).^2);
%GradB0 = @(x,y) -2*exp(-(7.5*x-5).^2-(7.5*y-5).^2)*[7.5*(7.5*x-5),7.5*(7.5*y-5)]...
%    +2*exp(-(7.5*x-2.5).^2-(7.5*y-2.5).^2)*[7.5*(7.5*x-2.5),7.5*(7.5*y-2.5)];

%Ground truth 3
B0 = @(x,y) exp(-(7.5*x-5.5).^2-(7.5*y-5.5).^2)...
    +.75*exp(-(5*x-1.25).^2-(7.5*y-5.5).^2)...
    +1.25*exp(-(7.5*x-5.5).^2-(5*y-1.25).^2)...
    +exp(-(7.5*x-2).^2-(7.5*y-2).^2);
GradB0 = @(x,y) -2*exp(-(7.5*x-5.5).^2-(7.5*y-5.5).^2)*[7.5*(7.5*x-5.5),7.5*(7.5*y-5.5)]...
    -2*.75*exp(-(5*x-1.25).^2-(7.5*y-5.5).^2)*[5*(5*x-1.25),7.5*(7.5*y-5.5)]...
    -2*1.25*exp(-(7.5*x-5.5).^2-(5*y-1.25).^2)*[7.5*(7.5*x-5.5),5*(5*y-1.25)]...
    -2*exp(-(7.5*x-2).^2-(7.5*y-2).^2)*[7.5*(7.5*x-2),7.5*(7.5*y-2)];

B0_bary=B0(barycenters(1,:),barycenters(2,:));
B0_norm = sqrt(sum(B0_bary.^2.*mesh_elements_area));
    % L^2 norm of B_0
B0_int = sum(B0_bary.*mesh_elements_area);

% Fix constant by setting integral to zero
B0 = @(x,y) B0(x,y) - B0_int;

% Evaluate B_0 at the mesh nodes and barycenters
B0_mesh=B0(mesh_nodes(1,:),mesh_nodes(2,:));
B0_bary=B0(barycenters(1,:),barycenters(2,:));

% Plot B_0
figure()
axes('FontSize', 20, 'NextPlot','add')
pdeplot(model,'XYData',real(B0_mesh),'ColorMap',jet)
axis equal
colorbar('Fontsize',20)
% pdeplot(model,'XYData',B0_mesh,'ZData',B0_mesh,'ColorMap',jet) 
    % 3D plot
title('True potential B_0','FontSize',20);
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
crameri vik


%%
% 2. Projection of B_0 onto the Fourier

% Specify truncation level
M=4;
param_num = 4*(M+1)^2;
    % number of basis functions, (M+1)^2 for each cos(2*pi*m*x)cos(2*pi*n*y), ...

% Specify the matrix of deterministic multiplicative coefficients
Kappa=4*ones(M+1,M+1);
Kappa(1,1)=1;
Kappa(1,2:M+1)=2;
Kappa(2:M+1,1)=2;

% Initialise matrices to store the Fourier coefficients
B0_cc_coeff=zeros(M+1,M+1);
B0_cs_coeff=zeros(M+1,M+1);
B0_sc_coeff=zeros(M+1,M+1);
B0_ss_coeff=zeros(M+1,M+1);

% Compute the Fourier coefficients
for m=0:M
    for n=0:M
        % Coefficients with respect to the CosCos functions
        emn_interp=scatteredInterpolant(mesh_nodes(1,:)',...
        mesh_nodes(2,:)',CosCos(m,n,mesh_nodes(1,:),mesh_nodes(2,:))');
        emn_bary=emn_interp(barycenters(1,:),barycenters(2,:));
        B0_cc_coeff(m+1,n+1)=sum(mesh_elements_area.*B0_bary.*emn_bary);

        % Coefficients with respect to the CosSin functions
        emn_interp=scatteredInterpolant(mesh_nodes(1,:)',...
        mesh_nodes(2,:)',CosSin(m,n,mesh_nodes(1,:),mesh_nodes(2,:))');
        emn_bary=emn_interp(barycenters(1,:),barycenters(2,:));
        B0_cs_coeff(m+1,n+1)=sum(mesh_elements_area.*B0_bary.*emn_bary);

        % Coefficients with respect to the SinCos functions
        emn_interp=scatteredInterpolant(mesh_nodes(1,:)',...
        mesh_nodes(2,:)',SinCos(m,n,mesh_nodes(1,:),mesh_nodes(2,:))');
        emn_bary=emn_interp(barycenters(1,:),barycenters(2,:));
        B0_sc_coeff(m+1,n+1)=sum(mesh_elements_area.*B0_bary.*emn_bary);

        % Coefficients with respect to the SinSin functions
        emn_interp=scatteredInterpolant(mesh_nodes(1,:)',...
        mesh_nodes(2,:)',SinSin(m,n,mesh_nodes(1,:),mesh_nodes(2,:))');
        emn_bary=emn_interp(barycenters(1,:),barycenters(2,:));
        B0_ss_coeff(m+1,n+1)=sum(mesh_elements_area.*B0_bary.*emn_bary);
    end
end

% Multipication by the deterministic coefficients
B0_cc_coeff=B0_cc_coeff.*Kappa;
B0_cs_coeff=B0_cs_coeff.*Kappa;
B0_sc_coeff=B0_sc_coeff.*Kappa;
B0_ss_coeff=B0_ss_coeff.*Kappa;

% Compute projection of B_0
B0_proj=zeros(1,mesh_nodes_num);
for m=0:M
    for n=0:M
        B0_proj = B0_proj...
            +B0_cc_coeff(m+1,n+1)*CosCos(m,n,mesh_nodes(1,:),mesh_nodes(2,:))...
            +B0_cs_coeff(m+1,n+1)*CosSin(m,n,mesh_nodes(1,:),mesh_nodes(2,:))...
            +B0_sc_coeff(m+1,n+1)*SinCos(m,n,mesh_nodes(1,:),mesh_nodes(2,:))...
            +B0_ss_coeff(m+1,n+1)*SinSin(m,n,mesh_nodes(1,:),mesh_nodes(2,:));
    end
end

% Plot B_0, its projection and the approximation error
figure()
subplot(1,2,1)
pdeplot(model,'XYData',real(B0_mesh),'ColorMap',jet)
title('True B_0','FontSize',20)
%clim([min(B0_mesh)+.5,max(B0_mesh)+.5])
colorbar('Fontsize',20)

subplot(1,2,2)
pdeplot(model,'XYData',real(B0_proj),'ColorMap',jet)
title('Projection of B_0','FontSize',20)
%clim([min(B0_mesh)+.5,max(B0_mesh)+.5])
colorbar('Fontsize',20)

% Interpolation of the Fourier projection at the barycenters
B0_proj_bary = griddata(mesh_nodes(1,:),mesh_nodes(2,:),B0_proj,...
    barycenters(1,:),barycenters(2,:));

% Approximate L^2 distance between B_0 and projection
approx_error = sqrt(sum((B0_bary-B0_proj_bary).^2.*mesh_elements_area));
disp(['L^2 approximation error via projection = ', num2str(approx_error)])
disp(['Relative error = ', num2str(approx_error/B0_norm)])

%%
% Repeated experiments to evaluate coverage for functionals

% Square functional
true_square_functional = sum((B0_bary.^2).*mesh_elements_area);
coverage_square = 0;
% q-power functional
q=4;
true_power_functional = sum((B0_bary.^4).*mesh_elements_area);
coverage_power = 0;
% Entrpy of invariant measure functional
mu0_bary=exp(B0_bary);
mu0_int = sum(mu0_bary.*mesh_elements_area);
mu0_bary = mu0_bary/mu0_int;
true_entropy_functional = sum(mu0_bary.*log(mu0_bary).*mesh_elements_area);
coverage_entropy = 0;

% Posterior draws to approximate one dimensional posteriors
post_draws_num = 1000;
post_draws_square_functional = zeros(1,post_draws_num);
post_draws_power_functional = zeros(1,post_draws_num);
post_draws_entropy_functional = zeros(1,post_draws_num);

% Prior covariance matrix
prior_cov=eye(4*(M+1)^2,4*(M+1)^2);
for i=1:(M+1)^2
        prior_cov(i,i)=1/(1+i);
        prior_cov((M+1)^2+i,(M+1)^2+i)=1/(1+i);
end

% Initialisation of diffusion paths
sigma=1;
sim_len = 100000;
    % number of simulated values
deltaT = .001;
    % small time interval between Euler-Maruyama iterates
T = sim_len *deltaT;
    % time horizon
X = zeros(2,sim_len); 
    % initialisation of the trajectory

% Number of repeated experiments
n_exp=50;
errors_square = zeros(1,n_exp);
errors_power = zeros(1,n_exp);
errors_entropy = zeros(1,n_exp);
length_square = zeros(1,n_exp);
length_power = zeros(1,n_exp);
length_entropy = zeros(1,n_exp);

% Initialisations
Sigma=ones((M+1)^2,(M+1)^2);
xderiv_cc=ones((M+1)^2,sim_len+1);
yderiv_cc=ones((M+1)^2,sim_len+1);
xderiv_cs=ones((M+1)^2,sim_len+1);
yderiv_cs=ones((M+1)^2,sim_len+1);
xderiv_sc=ones((M+1)^2,sim_len+1);
yderiv_sc=ones((M+1)^2,sim_len+1);
xderiv_ss=ones((M+1)^2,sim_len+1);
yderiv_ss=ones((M+1)^2,sim_len+1);
H=ones((M+1)^2,1);
post_mean_cc_matr=zeros((M+1),(M+1));
post_mean_cs_matr=zeros((M+1),(M+1));
post_mean_sc_matr=zeros((M+1),(M+1));
post_mean_ss_matr=zeros((M+1),(M+1));

rng(1)

for s=1:n_exp
    %disp(["Experiment n. ",num2str(s)])
    s

    % Sample trajectory
    x_0 = [rand();rand()];
    X(:,1) = x_0;
    X_prec=X(:,1);
    X_prec_proj=X_prec;
    for i=2:sim_len+1
        if X(1,i-1)>0
            X_prec_proj(1)=X(1,i-1)-fix(X(1,i-1));
        else
            X_prec_proj(1)=1+(X(1,i-1)-fix(X(1,i-1)));
        end

        if X(2,i-1)>0
            X_prec_proj(2)=X(2,i-1)-fix(X(2,i-1));
        else
            X_prec_proj(2)=1+(X(2,i-1)-fix(X(2,i-1)));
        end

        X(:,i) = X(:,i-1) - GradB0(X_prec_proj(1),X_prec_proj(2))'*deltaT...
            + sigma*mvnrnd(zeros(2,1),deltaT*eye(2))';
    end

    % Compute covariance matrix
    %Sigma=ones((M+1)^2,(M+1)^2);
        % initialises the matrix Sigma in the posterior covariance formula

    %xderiv_cc=ones((M+1)^2,sim_len+1);
    %yderiv_cc=ones((M+1)^2,sim_len+1);
    %xderiv_cs=ones((M+1)^2,sim_len+1);
    %yderiv_cs=ones((M+1)^2,sim_len+1);
    %xderiv_sc=ones((M+1)^2,sim_len+1);
    %yderiv_sc=ones((M+1)^2,sim_len+1);
    %xderiv_ss=ones((M+1)^2,sim_len+1);
    %yderiv_ss=ones((M+1)^2,sim_len+1);

    for m=0:M
        for n=0:M
        grad=GradCosCos(m,n,X(1,:),X(2,:));
        xderiv_cc((M+1)*m+(n+1),:)=grad(1,:);
        yderiv_cc((M+1)*m+(n+1),:)=grad(2,:);

        grad=GradCosSin(m,n,X(1,:),X(2,:));
        xderiv_cs((M+1)*m+(n+1),:)=grad(1,:);
        yderiv_cs((M+1)*m+(n+1),:)=grad(2,:);

        grad=GradSinCos(m,n,X(1,:),X(2,:));
        xderiv_sc((M+1)*m+(n+1),:)=grad(1,:);
        yderiv_sc((M+1)*m+(n+1),:)=grad(2,:);

        grad=GradSinSin(m,n,X(1,:),X(2,:));
        xderiv_ss((M+1)*m+(n+1),:)=grad(1,:);
        yderiv_ss((M+1)*m+(n+1),:)=grad(2,:);
        end
    end

    for i=1:(M+1)^2
        for j=1:(M+1)^2

        %cc-cc
        Sigma(i,j)=deltaT*sum(xderiv_cc(i,:).*xderiv_cc(j,:)...
            +yderiv_cc(i,:).*yderiv_cc(j,:));

        % cc-cs
        Sigma(i,(M+1)^2+j)=deltaT*sum(xderiv_cc(i,:).*xderiv_cs(j,:)...
            +yderiv_cc(i,:).*yderiv_cs(j,:));

        % cc-sc
        Sigma(i,2*(M+1)^2+j)=deltaT*sum(xderiv_cc(i,:).*xderiv_sc(j,:)...
            +yderiv_cc(i,:).*yderiv_sc(j,:));

        % cc-ss
        Sigma(i,3*(M+1)^2+j)=deltaT*sum(xderiv_cc(i,:).*xderiv_ss(j,:)...
            +yderiv_cc(i,:).*yderiv_ss(j,:));

        %%%%

        %cs-cc
        Sigma((M+1)^2+i,j)=deltaT*sum(xderiv_cs(i,:).*xderiv_cc(j,:)...
            +yderiv_cs(i,:).*yderiv_cc(j,:));

        % cs-cs
        Sigma((M+1)^2+i,(M+1)^2+j)=deltaT*sum(xderiv_cs(i,:).*xderiv_cs(j,:)...
            +yderiv_cs(i,:).*yderiv_cs(j,:));

        % cs-sc
        Sigma((M+1)^2+i,2*(M+1)^2+j)=deltaT*sum(xderiv_cs(i,:).*xderiv_sc(j,:)...
            +yderiv_cs(i,:).*yderiv_sc(j,:));

        % cs-ss
        Sigma((M+1)^2+i,3*(M+1)^2+j)=deltaT*sum(xderiv_cs(i,:).*xderiv_ss(j,:)...
            +yderiv_cs(i,:).*yderiv_ss(j,:));

        %%%

        %sc-cc
        Sigma(2*(M+1)^2+i,j)=deltaT*sum(xderiv_sc(i,:).*xderiv_cc(j,:)...
            +yderiv_sc(i,:).*yderiv_cc(j,:));

        % sc-cs
        Sigma(2*(M+1)^2+i,(M+1)^2+j)=deltaT*sum(xderiv_sc(i,:).*xderiv_cs(j,:)...
            +yderiv_sc(i,:).*yderiv_cs(j,:));

        % sc-sc
        Sigma(2*(M+1)^2+i,2*(M+1)^2+j)=deltaT*sum(xderiv_sc(i,:).*xderiv_sc(j,:)...
            +yderiv_sc(i,:).*yderiv_sc(j,:));

        % sc-ss
        Sigma(2*(M+1)^2+i,3*(M+1)^2+j)=deltaT*sum(xderiv_sc(i,:).*xderiv_ss(j,:)...
            +yderiv_sc(i,:).*yderiv_ss(j,:));

        %%%

        %ss-cc
        Sigma(3*(M+1)^2+i,j)=deltaT*sum(xderiv_ss(i,:).*xderiv_cc(j,:)...
            +yderiv_ss(i,:).*yderiv_cc(j,:));

        % ss-cs
        Sigma(3*(M+1)^2+i,(M+1)^2+j)=deltaT*sum(xderiv_ss(i,:).*xderiv_cs(j,:)...
            +yderiv_ss(i,:).*yderiv_cs(j,:));

        % ss-sc
        Sigma(3*(M+1)^2+i,2*(M+1)^2+j)=deltaT*sum(xderiv_ss(i,:).*xderiv_sc(j,:)...
            +yderiv_ss(i,:).*yderiv_sc(j,:));

        % ss-ss
        Sigma(3*(M+1)^2+i,3*(M+1)^2+j)=deltaT*sum(xderiv_ss(i,:).*xderiv_ss(j,:)...
            +yderiv_ss(i,:).*yderiv_ss(j,:));

        end
    end

    post_cov=inv(inv(prior_cov)+Sigma);
    post_cov(1,1)=0;

    % Compute posterior mean
    X_diff=X(:,2:sim_len)-X(:,1:sim_len-1);
    %H=ones((M+1)^2,1);
        % initialises vector H for computation of posterior mean

    for i=1:(M+1)^2
        H(i)=sum(xderiv_cc(i,1:sim_len-1).*X_diff(1,:))...
        +sum(yderiv_cc(i,1:sim_len-1).*X_diff(2,:));

        H((M+1)^2+i)=sum(xderiv_cs(i,1:sim_len-1).*X_diff(1,:))...
        +sum(yderiv_cs(i,1:sim_len-1).*X_diff(2,:));

        H(2*(M+1)^2+i)=sum(xderiv_sc(i,1:sim_len-1).*X_diff(1,:))...
        +sum(yderiv_sc(i,1:sim_len-1).*X_diff(2,:));

        H(3*(M+1)^2+i)=sum(xderiv_ss(i,1:sim_len-1).*X_diff(1,:))...
        +sum(yderiv_ss(i,1:sim_len-1).*X_diff(2,:));
    end

    post_mean=post_cov*H;

    % Posterior mean as a function
    %post_mean_cc_matr=zeros((M+1),(M+1));
    %post_mean_cs_matr=zeros((M+1),(M+1));
    %post_mean_sc_matr=zeros((M+1),(M+1));
    %post_mean_ss_matr=zeros((M+1),(M+1));
    %for m=0:M
    %    for n=0:M
    %        post_mean_cc_matr(m+1,n+1)=post_mean((M+1)*m+(n+1));
    %        post_mean_cs_matr(m+1,n+1)=post_mean((M+1)^2+(M+1)*m+(n+1));
    %        post_mean_sc_matr(m+1,n+1)=post_mean(2*(M+1)^2+(M+1)*m+(n+1));
    %        post_mean_ss_matr(m+1,n+1)=post_mean(3*(M+1)^2+(M+1)*m+(n+1));
    %    end
    %end
    %B0_mean_mesh=zeros(1,mesh_nodes_num);
    %for m=0:M
    %    for n=0:M
    %        B0_mean_mesh = B0_mean_mesh...
    %            +post_mean_cc_matr(m+1,n+1)*CosCos(m,n,mesh_nodes(1,:),mesh_nodes(2,:))...
    %            +post_mean_cs_matr(m+1,n+1)*CosSin(m,n,mesh_nodes(1,:),mesh_nodes(2,:))...
    %            +post_mean_sc_matr(m+1,n+1)*SinCos(m,n,mesh_nodes(1,:),mesh_nodes(2,:))...
    %            +post_mean_ss_matr(m+1,n+1)*SinSin(m,n,mesh_nodes(1,:),mesh_nodes(2,:));
    %    end
    %end
    %B0_mean_bary = griddata(mesh_nodes(1,:),mesh_nodes(2,:),B0_mean_mesh,...
    %barycenters(1,:),barycenters(2,:));
    %estim_error = sqrt(sum((B0_bary-sum(B0_bary.*mesh_elements_area)*ones(1,mesh_elements_num)-B0_mean_bary).^2.*mesh_elements_area));
    %disp(['Relative L^2 estimation error = ', num2str(estim_error/B0_norm)])

    % Posterior draws of functionals
    post_draws = mvnrnd(post_mean,post_cov,post_draws_num)';
    for j=1:post_draws_num
        post_draw=post_draws(:,j);
        % Exctract Fourier coefficients
        post_draw_cc_matr=zeros((M+1),(M+1));
        post_draw_cs_matr=zeros((M+1),(M+1));
        post_draw_sc_matr=zeros((M+1),(M+1));
        post_draw_ss_matr=zeros((M+1),(M+1));
        for m=0:M
            for n=0:M
                post_draw_cc_matr(m+1,n+1)=post_draw((M+1)*m+(n+1));
                post_draw_cs_matr(m+1,n+1)=post_draw((M+1)^2+(M+1)*m+(n+1));
                post_draw_sc_matr(m+1,n+1)=post_draw(2*(M+1)^2+(M+1)*m+(n+1));
                post_draw_ss_matr(m+1,n+1)=post_draw(3*(M+1)^2+(M+1)*m+(n+1));
            end
        end
        %B_draw_mesh=zeros(1,mesh_nodes_num);
        %for m=0:M
        %    for n=0:M
        %        B_draw_mesh = B_draw_mesh...
        %            +post_draw_cc_matr(m+1,n+1)*CosCos(m,n,mesh_nodes(1,:),mesh_nodes(2,:))...
        %            +post_draw_cs_matr(m+1,n+1)*CosSin(m,n,mesh_nodes(1,:),mesh_nodes(2,:))...
        %            +post_draw_sc_matr(m+1,n+1)*SinCos(m,n,mesh_nodes(1,:),mesh_nodes(2,:))...
        %            +post_draw_ss_matr(m+1,n+1)*SinSin(m,n,mesh_nodes(1,:),mesh_nodes(2,:));
        %    end
        %end
        %B_draw_bary = griddata(mesh_nodes(1,:),mesh_nodes(2,:),B_draw_mesh,...
        %barycenters(1,:),barycenters(2,:));
        B_draw_bary=zeros(1,mesh_elements_num);
        for m=0:M
            for n=0:M
                B_draw_bary = B_draw_bary...
                    +post_draw_cc_matr(m+1,n+1)*CosCos(m,n,barycenters(1,:),barycenters(2,:))...
                    +post_draw_cs_matr(m+1,n+1)*CosSin(m,n,barycenters(1,:),barycenters(2,:))...
                    +post_draw_sc_matr(m+1,n+1)*SinCos(m,n,barycenters(1,:),barycenters(2,:))...
                    +post_draw_ss_matr(m+1,n+1)*SinSin(m,n,barycenters(1,:),barycenters(2,:));
            end
        end
        % Posterior functionals
        post_draws_square_functional(j) = sum((B_draw_bary.^2).*mesh_elements_area);
        post_draws_power_functional(j) = sum((B_draw_bary.^q).*mesh_elements_area);
            % standard posterior
        %post_draws_functional(j) = sqrt(T)*(sum((B_draw_bary.^q).*mesh_elements_area)-true_power_functional);
            %scaled and centred posterior
        mu_draw_bary=exp(B_draw_bary);
        mu_draw_int = sum(mu_draw_bary.*mesh_elements_area);
        mu_draw_bary = mu_draw_bary/mu_draw_int;
        post_draws_entropy_functional(j) = sum(mu_draw_bary.*log(mu_draw_bary).*mesh_elements_area);
    end

    % Posterior mean estimation errors
    errors_square(s) = sqrt((true_square_functional-mean(post_draws_square_functional))^2);
    errors_power(s) = sqrt((true_power_functional-mean(post_draws_power_functional))^2);
    errors_entropy(s) = sqrt((true_entropy_functional-mean(post_draws_entropy_functional))^2);

    % Length of credible intervals
    length_square = quantile(post_draws_square_functional,0.975) - quantile(post_draws_square_functional,0.025);
    length_power = quantile(post_draws_power_functional,0.975) - quantile(post_draws_power_functional,0.025);
    length_entropy = quantile(post_draws_entropy_functional,0.975) - quantile(post_draws_entropy_functional,0.025);

    % Coverage checks
    if quantile(post_draws_square_functional,0.025) < true_square_functional
        if true_square_functional < quantile(post_draws_square_functional,0.975)
            coverage_square = coverage_square+1;
        end
    end
    if quantile(post_draws_power_functional,0.025) < true_power_functional
        if true_power_functional < quantile(post_draws_power_functional,0.975)
            coverage_power = coverage_power+1;
        end
    end
    if quantile(post_draws_entropy_functional,0.025) < true_entropy_functional
        if true_entropy_functional < quantile(post_draws_entropy_functional,0.975)
            coverage_entropy = coverage_entropy+1;
        end
    end
end

disp('###### COVERAGES ######')
coverage_square = coverage_square/n_exp;
disp(['Coverage square functional = ', num2str(coverage_square)])
coverage_power= coverage_power/n_exp;
disp(['Coverage power functional = ', num2str(coverage_power)])
coverage_entropy= coverage_entropy/n_exp;
disp(['Coverage entropy functional = ', num2str(coverage_entropy)])

disp('###### LENGTHS ######')
disp(['Length square functional = ', num2str(mean(length_square))])
disp(['Length power functional = ', num2str(mean(length_power))])
disp(['Length entropy functional = ', num2str(mean(length_entropy))])

disp('###### POSTERIOR MEAN ######')
disp(['RMSE square functional = ', num2str(mean(errors_square))])
disp(['STD DEV square functional = ', num2str(std(errors_square))])
disp(['RMSE power functional = ', num2str(mean(errors_power))])
disp(['STD DEV power functional = ', num2str(std(errors_power))])
disp(['RMSE entropy functional = ', num2str(mean(errors_entropy))])
disp(['STD DEV entropy functional = ', num2str(std(errors_entropy))])