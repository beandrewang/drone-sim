// an unit test for the kalman filter implementation
// an example of a velocity and position system

exec('kf.sci');

duration = 10;
dt = .1;

A = [1 dt; 0 1];
B = [dt ^ 2 / 2; dt];
C = [1 0];

u = 1.5;
X0 = [0; 0]
X = X0;
accel_noise_var = 0.05;
pos_noise_var = 10;
Ez = pos_noise_var ^ 2;
Ex = accel_noise_var ^ 2 * [dt ^ 4 / 4 dt ^ 3 / 2; dt ^ 3 / 2 dt ^ 2];
P = Ex;

Position = [];
Position_est = [];
vel = [];

subplot(211);
for t = 0 : dt : 10
    accel_noise = accel_noise_var * [(dt ^ 2 / 2 * rand(1, "normal")); dt * rand(1, "normal")];
    X = A * X + B * u + accel_noise;
    
    pos_noise = pos_noise_var * rand(1, "normal");
    y = C * X + pos_noise;
    Position = [Position; X(1)];
    Position_est = [Position_est, y];
    vel = [vel; X(2)];
    
    plot(0 : dt : t, Position, '-r.');
    plot(0 : dt : t, Position_est, '-k.');
end

Position_out = [];
Vel_out = [];
kf = kalman_filter_init(A, B, C, X0, P, u, Ex, Ez);
for i = 1 : length(Position)
    kf = kalman_filter_predict(kf);
    kf = kalman_filter_update(kf, Position_est(i));
    Position_out = [Position_out; kf.X(1)];
    Vel_out = [Vel_out; kf.X(2)];
end

subplot(212);
tt = 0 : dt : duration;
plot(tt, Position, '-r.');
plot(tt, Position_est, '-k.');
plot(tt, Position_out, '-g.');
