%% Masses
M1 = 1.55e-6;
M2 = 2.7e-6;
M3 = 25e-6;
M4 = 25e-6;
M5 = 1.78e-6;
M6 = 25.5e-6;

M = diag([M1, M2, M3, M4, M5, M6]);

%% Damping coefficients
C1 = 0.00007;
C2 = 0.5;
C3 = 1.74;
C4 = 0.122;
C5 = 0.359;
C6 = 0.00028;
C7 = 0.02;
C8 = 0.00004;
C9 = 0.1;
C10 = 0.1;

% Construct damping matrix C
C = zeros(6,6);

C(1,1) = C2;
C(1,2) = -C2;
C(2,1) = -C2;
C(2,2) = C1 + C2 + C3;
C(2,3) = -C3;
C(3,2) = -C3;
C(3,3) = C3 + C4 + C5;
C(3,4) = -C5;
C(4,3) = -C5;
C(4,4) = C5 + C6 + C7;
C(4,5) = -C6;
C(5,4) = -C6;
C(5,5) = C6 + C8 + C9;
C(5,6) = -C9;
C(6,5) = -C9;
C(6,6) = C9 + C10;

%% Stiffness coefficients
K1 = 1175;
K2 = 20001;
K3 = 94740;
K4 = 0;
K5 = 1000017; 
K6 = 167;
K7 = 0;
K8 = 623;

% Construct stiffness matrix K
K = zeros(6,6);

K(1,1) = K2;
K(1,2) = -K2;
K(2,1) = -K2;
K(2,2) = K1 + K2 + K3;
K(2,3) = -K3;
K(3,2) = -K3;
K(3,3) = K3 + K4 + K5;
K(3,4) = -K5;
K(4,3) = -K5;
K(4,4) = K5 + K6;
K(4,5) = -K6;
K(5,4) = -K6;
K(5,5) = K6 + K8;
K(5,6) = 0;
K(6,5) = 0;
K(6,6) = 0; % Assuming no stiffness associated with M6

%% Input matrix B
B = [0; 5.364e-5; 0; 0; 0; 0];

%% Output matrix H
H = [0 1 0 0 0 0; 
     0 0 0 0 1 0];

%% Loop over frequencies
% Define frequency range (in Hz)
f_range = logspace(log10(100), log10(10000), 1000); % 100 Hz to 10,000 Hz
omega = 2 * pi * f_range; % Convert to rad/s

% Initialize frequency response functions
G = zeros(2, length(omega));

% Loop over frequencies
for k = 1:length(omega)
    w = omega(k);
    sys_matrix = -M * w^2 + 1i * C * w + K;
    x = sys_matrix \ B;
    G(:, k) = H * x;
end

%% Convert magnitude to decibels
G_amp = abs(G);
G_amp = G_amp * 1e+6; % m to um
G_phase = rad2deg(angle(G));
% Subtract 360 degrees from phase angles greater than 0 degrees
G_phase(G_phase > 0) = G_phase(G_phase > 0) - 360;
% Plot the frequency response functions
figure;
subplot(2,1,1)
loglog(squeeze(f_range), squeeze(G_amp(1,:)),"b",squeeze(f_range), squeeze(G_amp(2,:)), "r");
xlabel('Frequency (Hz)');
ylabel('Magnitude (um)');
grid on;
subplot(2,1,2)
semilogx(squeeze(f_range), squeeze(G_phase(1,:)),"b",squeeze(f_range), squeeze(G_phase(2,:)), "r");
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
grid on;

%% Solve the generalized eigenvalue problem
[V, D] = eig(K, -M);

% The eigenvalues are on the diagonal of D
lambda = diag(D);

% Compute the natural frequencies (ωi) from λ = -ω^2
omega = sqrt(abs(lambda));

% Sort the natural frequencies and corresponding mode shapes
[omega_sorted, index] = sort(omega);
V_sorted = V(:, index);

% Display the results
disp('Natural Frequencies (rad/s) in Ascending Order:');
disp(omega_sorted);

disp('Mode Shapes (Eigenvectors) Corresponding to Each Frequency:');
disp(V_sorted);

%% Define the unperturbed positions of each mass
x_positions = [1, 4, 7, 10, 13, 16];
y_position = 1;  % Fixed y-position for all masses

% Plot each mode shape
figure;
for mode = 1:size(V_sorted, 2)
    % Create a new subplot for each mode shape
    subplot(2, 3, mode);  % Adjust to fit 6 mode shapes in a grid
    hold on;
    title(['Mode Shape ', num2str(mode)]);
    xlabel('Position');
    ylabel('Deflection');


    % Plot the deflected positions
    for i = 1:length(x_positions)
        deflection = V_sorted(i, mode);  % Scale the deflection
        rectangle('Position', [x_positions(i) + deflection, y_position, 1, 1], ...
                  'EdgeColor', 'red',LineWidth=1.5);
    end

    % Plot the initial unperturbed positions
    for i = 1:length(x_positions)
        rectangle('Position', [x_positions(i), y_position, 1, 1], ...
                  'EdgeColor', 'black', 'LineStyle', '--');
    end

    axis equal;
    xlim([0, 18]);  % Adjust limits to fit all masses
    ylim([-5, 5]);  % Adjust limits based on deflection
    hold off;
end
