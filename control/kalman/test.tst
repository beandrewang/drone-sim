// an unit test for the kalman filter implementation
// an example of a velocity and position system
// Pos = Pos + Vel * dt + adt^2 / 2;
// Vel = Vel + at;
// X = [Pos, Vel];
// y = Pos
// u = a
// X = [1, dt; 0 1]X + [dt^2 / 2; dt]u
// y = [1 0]X
// A = [1, dt; 0 1], B = [dt^2/2; dt], C = [1 0] 

exec('kf.sci');

duration = 10;
dt = .1;

A = [1 dt; 0 1];
B = [dt ^ 2 / 2; dt];
C = [1 0];

u = 1.5; // the accelerator
X0 = [0; 0]
X = X0;
accel_noise_var = 0.05; // the standard variance of the control noise
pos_noise_var = 10; // the standard varaiance of the measured noise
Ez = pos_noise_var ^ 2; // the variance of the control noise
Ex = accel_noise_var ^ 2 * [dt ^ 4 / 4 dt ^ 3 / 2; dt ^ 3 / 2 dt ^ 2]; // the covariance of the predict noise 
P = Ex; 

// generate the X with noise
Position = [];
Position_meas = []; // measured position
vel = [];
for t = 0 : dt : 10
    accel_noise = accel_noise_var * [(dt ^ 2 / 2 * rand(1, "normal")); dt * rand(1, "normal")];
    X = A * X + B * u + accel_noise;
    
    pos_noise = pos_noise_var * rand(1, "normal");
    y = C * X + pos_noise;
    Position = [Position; X(1)];
    Position_meas = [Position_meas, y];
    vel = [vel; X(2)];
end
subplot(211);
plot(0 : dt : t, Position, '-r.');
plot(0 : dt : t, Position_meas, '-k.');

// using kalman filter
Position_out = [];
Vel_out = [];
kf = kalman_filter_init(A, B, C, X0, P, u, Ex, Ez);
for i = 1 : length(Position)
    kf = kalman_filter_predict(kf);
    kf = kalman_filter_update(kf, Position_meas(i));
    Position_out = [Position_out; kf.X(1)];
    Vel_out = [Vel_out; kf.X(2)];
end
subplot(212);
tt = 0 : dt : duration;
plot(tt, Position, '-r.');
plot(tt, Position_meas, '-k.');
plot(tt, Position_out, '-g.');
