// This file implement a PID regulator
// o = Kp * e + Ki * E + Kd * (e' - e);

funcprot(0);

function p = pid_regulator_init(ref, feedback, dt, kp, ki, kd)
    // ref, the reference signal
    // feedback, the feedback signal
    // iout, the output of the integral
    // e, the last error
    // dt, the interval of pid regulate
    // kp, the proportion coefficient
    // ki, the intergration coefficient
    // kd, the derivative coefficient
    p = struct('ref', ref, 'feedback', feedback, 'iout', 0, 'e', 0, 'dt', ...
        dt, 'kp', kp, 'ki', ki, 'kd', kd, 'output', 0);
endfunction

function p = pid_regulator_update(p, feedback)
    // p, the pid structure
    e = (p.ref - p.feedback);
    p.feedback = feedback;
    p.iout = p.iout + e * p.dt;
    if 0 == p.e then
        p.output = p.kp * e + p.ki * p.iout;
    else
        p.output = p.kp * e + p.ki * p.iout + p.kd * (e - p.e) / p.dt;
    end
    
    p.e = e;
endfunction
