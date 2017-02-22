// implementation of svpwm

function [_alpha, _beta] = forward_clark_transform(a, b, c)
    // this is a implementation of power-invariant 
    // a, b, c, the 3 phase voltages
    // _alpha, the alpha component in alpha-beta-gamma space
    // _beta, the beta component in alpha-beta-gamma space
    _alpha = a;
    _beta  = (b - c) * (1 / sqrt(3));
endfunction

function [a, b, c] = reverse_clack_transform(_alpha, _beta)
    // a, b, c, the 3 phase voltages
    // _alpha, the alpha component in alpha-beta space
    // _beta, the beta component in alpha-beta space
    a = _alpha;
    b = -_alpha / 2 + sqrt(3) * _beta / 2;
    c = -_alpha / 2 - sqrt(3) * _beta / 2;
endfunction

function [d, q] = forward_park_transform(_alpha, _beta, _theta)
    // _alpha, the alpha component in alpha-beta space
    // _beta, the beta component in alpha-beta space
    // _theta, the angle that the current motor rotates
    // d, the d component of the motor, torque
    // q, the q component of the motor, flux
    
    d = _alpha * cos(_theta) + _beta * sin(_theta);
    q = -_alpha * sin(_theta);
endfunction

function [_alpha, _beta] = reverse_park_transform(d, q, _theta)
    // _alpha, the alpha component in alpha-beta space
    // _beta, the beta component in alpha-beta space
    // _theta, the angle that the current motor rotates
    // d, the d component of the motor
    // q, the q component of the motor
    
    _alpha = d * cos(_theta) - q * sin(_theta);
    _beta  = d * sin(_theta) + q * cos(_theta);
endfunction

function [a, b, c, t0, t1, t2, n, clock_freq] = svpwm_period_generator(_alpha, _beta, T, Vdc)
    // _alpha, the alpha component in alpha-beta phase
    // _beta, the beta component in alpha-beta phase
    // T, the switching perioid, usually, this is the time of the loop.
    // Vdc, the DC source voltage
    // t0, the zero voltage perioid
    // t1, the Vx voltage perioid
    // t2, the Vy voltage perioid
    // n, the current sector
    // a, the generated a control
    // b, the generated b control
    // c, the generated c control
    // clock_freq, the pwm source clock frequency
    
    m = sqrt(_alpha ^ 2 + _beta ^ 2) / (2 / 3 * Vdc); // should < 0.907. under modulation
    theta = atan(_beta, _alpha);
    
    if theta >= 0 & theta < %pi / 3  then
        // sector 1
        n = 1;
    elseif theta >= %pi / 3 & theta < %pi * 2 / 3 then
        // sector 2
        n = 2;
    elseif theta >= %pi * 2 / 3 & theta < %pi then
        // sector 3
        n = 3;
    elseif theta >= -%pi & theta < -%pi * 2 / 3 then
        // sector 4
        n = 4;
    elseif theta >= -%pi * 2 / 3 & theta < -%pi / 3 then
        // sector 5
        n = 5;
    elseif theta >= -%pi / 3 & theta < 0 then
        // sector 6
        n = 6;
    end
    
    t1 = T * m * sin(n * %pi / 3 - theta);
    t2 = T * m * sin(theta - (n - 1) * %pi / 3);
    t0 = T - t1 - t2;
    
    clock_freq = T / 100;
    scale = 1 / clock_freq;
    T = T * scale;
    T0 = t0 * scale;
    T1 = t1 * scale;
    T2 = t2 * scale;
    a = zeros(1, T + 1);
    b = zeros(1, T + 1);
    c = zeros(1, T + 1);
    
    if 1 == n then
        // sector 1
        a(1 + T0 / 4 : T - T0 / 4) = 1;
        b(1 + (T0 / 4 + T1 / 2) : T - (T0 / 4 + T1 / 2)) = 1;
        c(1 + (T0 / 4 + T1 / 2 + T2 / 2) : T - (T0 / 4 + T1 / 2 + T2 / 2)) = 1;
    elseif 2 == n then
        // sector 2
        a(1 + (T0 / 4 + T2 / 2) : T - (T0 / 4 + T2 / 2)) = 1;
        b(1 + T0 / 4 : T - T0 / 4) = 1;
        c(1 + (T0 / 4 + T1 / 2 + T2 / 2) : T - (T0 / 4 + T1 / 2 + T2 / 2)) = 1;
    elseif 3 == n then
        // sector 3
        a(1 + (T0 / 4 + T1 / 2 + T2 / 2) : T - (T0 / 4 + T1 / 2 + T2 / 2)) = 1;
        b(1 + T0 / 4 : T - T0 / 4) = 1;
        c(1 + (T0 / 4 + T1 / 2) : T - (T0 / 4 + T1 / 2)) = 1;
    elseif 4 == n then
        // sector 4
        a(1 + (T0 / 4 + T1 / 2 + T2 / 2) : T - (T0 / 4 + T1 / 2 + T2 / 2)) = 1;
        b(1 + (T0 / 4 + T2 / 2) : T - (T0 / 4 + T2 / 2)) = 1;
        c(1 + T0 / 4 : T - T0 / 4) = 1;
    elseif 5 == n then
        // sector 5
        a(1 + (T0 /4 + T1 / 2) : T - (T0 / 4 + T1 / 2)) = 1;
        b(1 + (T0 / 4 + T1 / 2 + T2 / 2) : T - (T0 / 4 + T1 / 2 + T2 / 2)) = 1;
        c(1 + T0 / 4 : T - T0 / 4) = 1;
    elseif 6 == n then
        // sector 6
        a(1 + T0 / 4 : T - T0 / 4) = 1;
        b(1 + (T0 / 4 + T1 / 2 + T2 / 2) : T - (T0 / 4 + T1 / 2 + T2 / 2)) = 1;
        c(1 + (T0 /4 + T2 / 2) : T - (T0 / 4 + T2 / 2)) = 1;
    end
endfunction

function [Sa, Sb, Sc] = plot_svpwm(t0, t1, t2, n, Vdc, start_time)
    // generate a prioid SVPWM curve
    // t0, the zero voltage perioid
    // t1, the Vx voltage perioid
    // t2, the Vy voltage perioid
    // n, current sector
    // Vdc, the DC source voltage
    // start_time, the start location you want to continue to plot
    
    clock_freq = (t0 + t1 + t2) / 1000;
    scale = 1 / clock_freq;
    T = ceil((t0 + t1 + t2) * scale);
    a = zeros(1, T + 1);
    b = zeros(1, T + 1);
    c = zeros(1, T + 1);
    A = zeros(1, T + 1);
    B = zeros(1, T + 1);
    C = zeros(1, T + 1);
    N = ones(1, T + 1);
    
    T0 = t0 * scale;
    T1 = t1 * scale;
    T2 = t2 * scale;
    
    if 1 == n then
        // sector 1
        a(1 + T0 / 4 : T - T0 / 4) = 1;
        b(1 + (T0 / 4 + T1 / 2) : T - (T0 / 4 + T1 / 2)) = 1;
        c(1 + (T0 / 4 + T1 / 2 + T2 / 2) : T - (T0 / 4 + T1 / 2 + T2 / 2)) = 1;
        
        //A(1 : T0 / 4) = 0;
        A(T0 / 4 + 1 : (T0 / 4 + T1 / 2)) = 2 / 3 * Vdc; 
        A(T0 / 4 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = 1 / 3 * Vdc;
        //A(T0/ / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2) = 1 / 3 * Vdc;
        A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = 2 / 3 * Vdc;
       // A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        //B(1 : T0 / 4) = 0;
        B(T0 / 4 + 1 : (T0 / 4 + T1 / 2)) = -1 / 3 * Vdc;
        B(T0 / 4 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = 1 / 3 * Vdc;
        //B(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2) = 1 / 3 * Vdc;
        B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = -1 / 3 * Vdc;
       // B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        //C(1 : T0 / 4) = 0;
        C(T0 / 4 + 1 : (T0 / 4 + T1 / 2)) = -1 / 3 * Vdc;
        C(T0 / 4 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = -2 / 3 * Vdc;
        //C(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2) = -2 / 3 * Vdc;
        C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = -1 / 3 * Vdc;
        // C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
       
       Sa = (T - T0 / 2) / T;
       Sb = (T - T0 / 2 - T1) / T;
       Sc = (T0 / 2) / T;
       
       SA = (T1 / T * 2 / 3 + T2 / T * 1 / 3) * Vdc;
       SB = (-T1 / T * 1 / 3 + T2 / T * 1 / 3) * Vdc;
       SC = (-T1 / T * 1 / 3 - T2 / T * 2 / 3) * Vdc;
       
       SAB = T1 / T * Vdc;
       SBC = T2 / T * Vdc;
       SCA = -(T1 + T2) / T * Vdc;
        
    elseif 2 == n then
        // sector 2
        a(1 + (T0 / 4 + T2 / 2) : T - (T0 / 4 + T2 / 2)) = 1;
        b(1 + T0 / 4 : T - T0 / 4) = 1;
        c(1 + (T0 / 4 + T1 / 2 + T2 / 2) : T - (T0 / 4 + T1 / 2 + T2 / 2)) = 1;
        
        //A(1 : T0 / 4) = 0;
        A(T0 / 4 + 1 : (T0 / 4 + T2 / 2)) = -1 / 3 * Vdc;
        A(T0 / 4 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = 1 / 3 * Vdc;
        //A(T0/ / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2) = 1 / 3 * Vdc;
        A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = -1 / 3 * Vdc;
       // A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
       
        //B(1 : T0 / 4) = 0;
        B(T0 / 4 + 1 : (T0 / 4 + T2 / 2)) = 2 / 3 * Vdc;
        B(T0 / 4 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = 1 / 3 * Vdc;
        //B(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2) = 1 / 3 * Vdc;
        B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = 2 / 3 * Vdc;
        //B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        //C(1 : T0 / 4) = 0;
        C(T0 / 4 + 1 : (T0 / 4 + T2 / 2)) = -1 / 3 * Vdc;
        C(T0 / 4 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = -2 / 3 * Vdc;
        //C(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2) = -2 / 3 * Vdc;
        C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = -1 / 3 * Vdc;
        //C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        Sa = (T - T0 / 2 - T2) / T;
        Sb = (T - T0 / 2) / T;
        Sc = T0 / 2 / T;
        
        SA = (-T2 / T * 1 / 3 + T1 / T * 1 / 3) * Vdc;
        SB = (T2 / T * 2 / 3 + T1 / T * 1 / 3) * Vdc;
        SC = (-T2 / T * 1 / 3 - T1 / T * 2 / 3) * Vdc;
        
        SAB = -T2 / T * Vdc;
        SBC = (T1 + T2) / T * Vdc;
        SCA = -T1 / T * Vdc;
    elseif 3 == n then
        // sector 3
        a(1 + (T0 / 4 + T1 / 2 + T2 / 2) : T - (T0 / 4 + T1 / 2 + T2 / 2)) = 1;
        b(1 + T0 / 4 : T - T0 / 4) = 1;
        c(1 + (T0 / 4 + T1 / 2) : T - (T0 / 4 + T1 / 2)) = 1;
        
        //A(1 : T0 / 4) = 0;
        A(T0 / 4 + 1 : (T0 / 4 + T1 / 2)) = -1 / 3 * Vdc;
        A(T0 / 4 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = -2 / 3 * Vdc;
        //A(T0/ / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2) = -2 / 3 * Vdc;
        A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = -1 / 3 * Vdc;
       // A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        //B(1 : T0 / 4) = 0;
        B(T0 / 4 + 1 : (T0 / 4 + T1 / 2)) = 2 / 3 * Vdc;
        B(T0 / 4 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = 1 / 3 * Vdc;
        //B(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2) = 1 / 3 * Vdc;
        B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = 2 / 3 * Vdc;
        //B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        //C(1 : T0 / 4) = 0;
        C(T0 / 4 + 1 : (T0 / 4 + T1 / 2)) = -1 / 3 * Vdc;
        C(T0 / 4 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = 1 / 3 * Vdc;
        //C(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2) = 1 / 3 * Vdc;
        C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = -1 / 3 * Vdc;
        //C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        Sa = T0 / 2 / T;
        Sb = (T - T0 / 2) / T;
        Sc = (T - T0 / 2 - T1) / T;
        
        SA = (-T1 / T * 1 / 3 - T2 / T * 2 / 3) * Vdc;
        SB = (T1 / T * 2 / 3 + T2 / T * 1 / 3) * Vdc;
        SC = (-T1 / T * 1 / 3 + T2 / T * 1 / 3) * Vdc;
        
        SAB = -(T1 + T2) / T * Vdc;
        SBC = T1 / T * Vdc;
        SCA = T2/ T * Vdc;
    elseif 4 == n then
        // sector 4
        a(1 + (T0 / 4 + T1 / 2 + T2 / 2) : T - (T0 / 4 + T1 / 2 + T2 / 2)) = 1;
        b(1 + (T0 / 4 + T2 / 2) : T - (T0 / 4 + T2 / 2)) = 1;
        c(1 + T0 / 4 : T - T0 / 4) = 1;
        
        //A(1 : T0 / 4) = 0;
        A(T0 / 4 + 1 : (T0 / 4 + T2 / 2)) = -1 / 3 * Vdc;
        A(T0 / 4 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = -2 / 3 * Vdc;
        //A(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2) = -2 / 3 * Vdc;
        A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = -1 / 3 * Vdc;
        //A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        //B(1 : T0 / 4) = 0;
        B(T0 / 4 + 1 : (T0 / 4 + T2 / 2)) = -1 / 3 * Vdc;
        B(T0 / 4 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = 1 / 3 * Vdc;
        //B(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2) = 1 / 3 * Vdc;
        B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = -1 / 3 * Vdc;
        //B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        //C(1 : T0 / 4) = 0;
        C(T0 / 4 + 1 : (T0 / 4 + T2 / 2)) = 2 / 3 * Vdc;
        C(T0 / 4 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = 1 / 3 * Vdc;
        //C(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2) = 1 / 3 * Vdc;
        C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = 2 / 3 * Vdc;
        //C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        Sa = T0 / 2 / T;
        Sc = (T - T0 / 2) / T;
        Sb = (T - T0 / 2 - T2) / T;
        
        SA = (-T2 / T * 1 / 3 - T1 / T * 2 / 3) * Vdc;
        SB = (-T2 / T * 1 / 3 + T1 / T * 1 / 3) * Vdc;
        SC = (T2 / T * 2 / 3 + T1 / T * 1 / 3) * Vdc;
        
        SAB = -(T1) / T * Vdc;
        SBC = -T2 / T * Vdc;
        SCA = (T1 + T2)/ T * Vdc;
    elseif 5 == n then
        // sector 5
        a(1 + (T0 /4 + T1 / 2) : T - (T0 / 4 + T1 / 2)) = 1;
        b(1 + (T0 / 4 + T1 / 2 + T2 / 2) : T - (T0 / 4 + T1 / 2 + T2 / 2)) = 1;
        c(1 + T0 / 4 : T - T0 / 4) = 1;
        
        //A(1 : T0 / 4) = 0;
        A(T0 / 4 + 1 : (T0 / 4 + T1 / 2)) = -1 / 3 * Vdc;
        A(T0 / 4 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = 1 / 3 * Vdc;
        //A(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2) = 1 / 3 * Vdc;
        A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = -1 / 3 * Vdc;
        //A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        //B(1 : T0 / 4) = 0;
        B(T0 / 4 + 1 : (T0 / 4 + T1 / 2)) = -1 / 3 * Vdc;
        B(T0 / 4 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = -2 / 3 * Vdc;
        B(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2) = -2 / 3 * Vdc;
        B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = -1 / 3 * Vdc;
        //B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        //C(1 : T0 / 4) = 0;
        C(T0 / 4 + 1 : (T0 / 4 + T1 / 2)) = 2 / 3 * Vdc;
        C(T0 / 4 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = 1 / 3 * Vdc;
        //C(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2) = 1 / 3 * Vdc;
        C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = 2 / 3 * Vdc;
        //C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        Sb = T0 / 2 / T;
        Sc = (T - T0 / 2) / T;
        Sa = (T - T0 / 2 - T1) / T;
        
        SA = (-T1 / T * 1 / 3 + T2 / T * 1 / 3) * Vdc;
        SB = (-T1 / T * 1 / 3 - T2 / T * 2 / 3) * Vdc;
        SC = (T1 / T * 2 / 3 + T2 / T * 1 / 3) * Vdc;
        
        SAB = (T2) / T * Vdc;
        SBC = -(T1 + T2) / T * Vdc;
        SCA = T1/ T * Vdc;
    elseif 6 == n then
        // sector 6
        a(1 + T0 / 4 : T - T0 / 4) = 1;
        b(1 + (T0 / 4 + T1 / 2 + T2 / 2) : T - (T0 / 4 + T1 / 2 + T2 / 2)) = 1;
        c(1 + (T0 /4 + T2 / 2) : T - (T0 / 4 + T2 / 2)) = 1;
        
        //A(1 : T0 / 4) = 0;
        A(T0 / 4 + 1 : (T0 / 4 + T2 / 2)) = 2 / 3 * Vdc;
        A(T0 / 4 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = 1 / 3 * Vdc;
        //A(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2) = 1 / 3 * Vdc;
        A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = 2 / 3 * Vdc;
        //A(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        //B(1 : T0 / 4) = 0;
        B(T0 / 4 + 1 : (T0 / 4 + T2 / 2)) = -1 / 3 * Vdc;
        B(T0 / 4 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = -2 / 3 * Vdc;
        //B(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2) = -2 / 3 * Vdc;
        B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = -1 / 3 * Vdc;
        //B(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        //C(1 : T0 / 4) = 0;
        C(T0 / 4 + 1 : (T0 / 4 + T2 / 2)) = -1 / 3 * Vdc;
        C(T0 / 4 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = 1 / 3 * Vdc;
        //C(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2) = 1 / 3 * Vdc;
        C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = -1 / 3 * Vdc;
        //C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        Sb = T0 / 2 / T;
        Sa = (T - T0 / 2) / T;
        Sc = (T - T0 / 2 - T2) / T;
        
        SA = (T2 / T * 2 / 3 + T1 / T * 1 / 3) * Vdc;
        SB = (-T2 / T * 1 / 3 - T1 / T * 2 / 3) * Vdc;
        SC = (-T2 / T * 1 / 3 + T1 / T * 1 / 3) * Vdc;
        
        SAB = (T1 + T2) / T * Vdc;
        SBC = -T1 / T * Vdc;
        SCA = -T2/ T * Vdc;
    end
        
    //disp(start_time : 1 / scale : start_time + T / scale);
    f1 = scf(1);
    subplot(211);
    plot(start_time : 1 / scale : start_time + T / scale, a + 2, 'r');
    plot(start_time : 1 / scale : start_time + T / scale, c - 2, 'b');
    plot(start_time : 1 / scale : start_time + T / scale, b, 'g');
    xtitle('svpwm curve', 'time', 'value');
    legend('a', 'b', 'c');
    
    subplot(212);
    plot(start_time : 1 / scale : start_time + T / scale, A, 'r');
    plot(start_time : 1 / scale : start_time + T / scale, B, 'g');
    plot(start_time : 1 / scale : start_time + T / scale, C, 'b');
    xtitle('source of inverter', 'time', 'value');
    legend('a', 'b', 'c');
    //plot(start_time : 1 / scale : start_time + T / scale, N * n, '*c');
    
    f2 = scf(2);
    subplot(311)
    ppSa = ones(1, T + 1) * Sa;
    plot(start_time : 1 / scale : start_time + T / scale, ppSa, 'r');
    ppSb = ones(1, T + 1) * Sb;
    plot(start_time : 1 / scale : start_time + T / scale, ppSb, 'g');
    ppSc = ones(1, T + 1) * Sc;
    plot(start_time : 1 / scale : start_time + T / scale, ppSc, 'b');
    xtitle('svpwm equivalent voltage', 'time', 'value');
    legend('a', 'b', 'c');
    
    subplot(312);
    ppSa = ones(1, T + 1) * SA;
    plot(start_time : 1 / scale : start_time + T / scale, ppSa, 'r');
    ppSb = ones(1, T + 1) * SB;
    plot(start_time : 1 / scale : start_time + T / scale, ppSb, 'g');
    ppSc = ones(1, T + 1) * SC;
    plot(start_time : 1 / scale : start_time + T / scale, ppSc, 'b');
    xtitle('inverter output - phase voltage', 'time', 'value');
    legend('phase a', 'phase b', 'phase c');
    
    subplot(313);
    ppSa = ones(1, T + 1) * SAB;
    plot(start_time : 1 / scale : start_time + T / scale, ppSa, 'r');
    ppSb = ones(1, T + 1) * SBC;
    plot(start_time : 1 / scale : start_time + T / scale, ppSb, 'g');
    ppSc = ones(1, T + 1) * SCA;
    plot(start_time : 1 / scale : start_time + T / scale, ppSc, 'b');
    xtitle('inverter output - line voltage', 'time', 'value');
    legend('ab', 'bc', 'ca');
endfunction


