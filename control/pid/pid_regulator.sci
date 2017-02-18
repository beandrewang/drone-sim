// This file implement a PID regulator
// o = Kp * e + Ki * E + Kd * (e' - e);

// global variants

function [o, E, e_prev] = pid_regulator(r, f, E, e_prev, dt, Kp, Ki, Kd)
    // r, the reference signal
    // f, the feedback signal
    // o, the output of the PID regulator
    // E, the output of the intergrator
    // e_prev, the previouse error signal
    // Kp, the proportion coefficient 
    // Ki, the intergration coefficient
    // Kd, the derivative coefficient
    
    e = r - f; // current error signal
    E = E + e * dt; // 
    if e_prev == 0 then
        o = Kp * e + Ki * E;
    else
        o = Kp * e + Ki * E + Kd * (e- e_prev) * dt;
    end
    e_prev = e;
    //mprintf('%d, %d, %d', e, E, o);
    //pause
endfunction
