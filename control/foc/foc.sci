// implementation of FOC - field oriented control
exec('../pid/pid_regulator.sci');
exec('../svpwm/svpwm_generator.sci');
exec('../inverter/inverter.sci');
exec('../motor/motor.sci');

function [Ud, E, e_prev] = d_pi(Id, Id_back, E, e_prev, dt, Kp, Ki)
    // Id, the expected Id component, also means the flux
    // Id_back, the feedback Id component from motor
    [Ud, E, e_prev] = pid_regulator(Id, Id_back, E, e_prev, dt, Kp, Ki, 0);
endfunction

function [Uq, E, e_prev] = q_pi(Iq, Iq_back, E, e_prev, dt, Kp, Ki)
    // Iq, the expected Iq component, also means the flux
    // Iq_back, the feedback Iq component from motor
    [Uq, E, e_prev] = pid_regulator(Iq, Iq_back, E, e_prev, dt, Kp, Ki, 0);
endfunction

function [Id_back, Iq_back, theta, Ed, ed_prev, Eq, eq_prev] = foc(Id, Iq, Vdc, R, theta, Id_back, Iq_back, Ed, Eq, ed, eq)
    // this is a simple implementation for a FOC, there are two PID regulator 
    // in the design, one for the Id component and the other for the Iq component. 
    // Id, the expected Id component, also means the flux
    // Iq, the expected Iq component, also means the torque
    // Id_back, the feedback Id component from motor
    // Iq_back, the feedback Iq component from motor
    // theta, the current 
    //
    
    Kdp = 1;
    Kdi = 0;
    Kqp = 1;
    Kqi = 0;
    T = 0.2;
    dt = 0.1;
   
    [Vd, Ed, ed_prev] = d_pi(Id, Id_back, Ed, ed, dt, Kdp, Kdi);
    [Vq, Eq, eq_prev] = q_pi(Iq, Iq_back, Eq, eq, dt, Kqp, Kqi);
    
    [Valpha, Vbeta] = reverse_park_transform(Vd, Vq, theta);
    [a, b, c, t0, t1, t2, n, clock_freq, m, theta] = svpwm_period_generator(Valpha, Vbeta, T, Vdc);
    [A, B, C] = l2p3_inverter(a, b, c, t0, t1, t2, n, Vdc);
    [Ia, Ib, Ic, M, theta] = pmsm(A, B, C, R);
    [Ialpha, Ibeta] = forward_clark_transform(Ia, Ib, Ic);
    [Id_back, Iq_back] = forward_park_transform(Ialpha, Ibeta, theta);
    
    if ed == 0 & eq == 0 then
        break;
    end
endfunction
