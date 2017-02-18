// This is a test file for pid

exec('pid_regulator.sci', -1);

ref = 100
dt = 0.1
out = zeros(100 / dt);
t = 1 : dt : 100;
m = 2;
v = 0;
e_prev = 0;
E = 0;

function v = velocity_controller (v, a, dt)
    // v, the current velocity
    // a, the accelarater 
    // dt, the refresh time
    v = v + a * dt;
endfunction

i = 1;

for j = t
    // try to change the Kp and Ki, 
    // Kp from 1 to 10
    // Ki from 0 to 2

    [u, E, e_prev] = pid_regulator(100, out(i), E, e_prev, dt, 1, 0, 0);
    out(i + 1) = velocity_controller(out(i), u);
    i = i + 1;
end

plot(t, out(1 : $ - 1));

