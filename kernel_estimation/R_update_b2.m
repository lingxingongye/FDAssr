function   R    =   R_update_b2(MSI_BS, HSI_2D,R,mu)
for i = 1:size(MSI_BS,1)
    y = MSI_BS(i,:)';
    x = HSI_2D';
    [m,n]=size(x);

% B-spline basis setup
k = 3; % B-spline degree (cubic B-spline)
num_knots = n+k; % Set number of knots equal to the number of samples to potentially overfit
knots = augknt(linspace(0, 1, num_knots - k + 2), k); % Define the knots for the B-spline

% Create B-spline basis matrix for phi
phi_basis = zeros(num_knots + k - 2, n); % Basis matrix for phi
for i = 1:n
     % Assume each dimension in phi is parameterized by a different B-spline
   sp = spmak(knots, (1:num_knots + k - 2 == i)); % Create B-spline basis function for each dimension
    phi_basis(:, i) = fnval(sp, linspace(0, 1, num_knots + k - 2)); % Evaluate the B-spline at appropriate points
end

% Form the linear model to estimate coefficients of B-spline representation
H = x * phi_basis'; % H is the design matrix for the regression problem

% Estimate the B-spline coefficients for phi
b_coeff = (H' * H) \ (H' * y); % Solve normal equations to estimate coefficients
% Reconstruct estimated phi using B-spline coefficients
phi_estimated = phi_basis * b_coeff;

R(i,:) = phi_estimated;
end
end

% ======================
% Function to compute B-spline basis matrix
function basis_values = bspline_basismatrix(degree, knots, x)
    % Compute the B-spline basis values for a given degree, knots, and input value x
    % degree: Degree of the B-spline
    % knots: Knot vector
    % x: Evaluation point
    basis_values = spcol(knots, degree + 1, x); % This uses the spline toolbox's spcol function
end
