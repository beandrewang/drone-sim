// an unit test for the kalman filter implementation
// an example of a velocity and position system

exec('kf.sci');

Pos = 0; // the initial position of the system
Vel = 10; // the mean velocity of the system
Pos_sigma = 0; // the invariant of the position
Vel_sigma = 1; // the invariant of the velocity
dt = 1; // the sample interval
// the system state equation is
// Pos = Pos_prev + dt * Vel_prev
// Vel = Vel_prev
// [Pos; Vel] = [1 dt; 0 1] * [Pos_prev; Vel_prev] 
// dX = [1 dt; 0 1]X
// 
A = [1 dt; 0 1];
B = 0;
X0 = [Pos; Vel];
P0 = [Pos_sigma ^ 2 Pos_sigma * Vel_sigma; Pos_sigma * Vel_sigma Vel_sigma ^ 2]; // the intial state coinvariant
Q = [1, 2; 2, 4]; // the predict coinvariant

mPos_sigma = 20;
mVel_sigma = 5;
R = [mPos_sigma ^ 2, mPos_sigma * mVel_sigma; mPos_sigma * mVel_sigma, mVel_sigma ^ 2];
H = [1 0; 0 1];  // the same unit

Vel_array = grand(1, 100, 'nor', Vel, mVel_sigma);
Pos_array = Vel .* dt .* [1 : 100] + grand(1, 100, 'nor', 0, mPos_sigma);

subplot(211);

subplot(212);


sPos_array = zeros(1, 100);
sVel_array = zeros(1, 100);
kf = kalman_filter_init(A, B, X0, P0, 0, Q, R, H);

for i = 1 : 100
    kf = kalman_filter_predict(kf);
    kf = kalman_filter_update(kf, [Pos_array(i); Vel_array(i)]);
    sPos_array(i) = kf.X(1);
    sVel_array(i) = kf.X(2);
end

subplot(211);
plot([1 : 100], Pos_array);
plot([1 : 100], sPos_array, 'r');
xtitle('position', 't', 'p');
legend('raw', 'smoothed');
subplot(212);
plot([1 : 100], Vel_array);
plot([1 : 100], sVel_array, 'r');
xtitle('velocity', 't', 'v');
legend('raw', 'smoothed');

