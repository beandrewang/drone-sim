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

function [t0, t1, t2, n] = svpwm_period_generator(_alpha, _beta, T, Vdc)
    // _alpha, the alpha component in alpha-beta phase
    // _beta, the beta component in alpha-beta phase
    // T, the switching perioid, usually, this is the time of the loop.
    // Vdc, the DC source voltage
    // t0, the zero voltage perioid
    // t1, the Vx voltage perioid
    // t2, the Vy voltage perioid
    // n, the current sector
    
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
endfunction

function [Sa, Sb, Sc] = plot_svpwm(t0, t1, t2, n, Vdc, start_time)
    // generate a prioid SVPWM curve
    // t0, the zero voltage perioid
    // t1, the Vx voltage perioid
    // t2, the Vy voltage perioid
    // n, current sector
    // Vdc, the DC source voltage
    // start_time, the start location you want to continue to plot
    
    scale = 100;
    
    T = (t0 + t1 + t2) * scale;
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
       
       Sa = (T - T0 / 2) / T * Vdc;
       Sb = (T - T0 / 2 - T1) / T * Vdc;
       Sc = (T0 / 2) / T * Vdc;s
        
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
        
        Sa = (T - T0 / 2 - T2) / T * Vdc;
        Sb = (T - T0 / 2) / T * Vdc;
        Sc = T0 / 2 / T * Vdc;
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
        
        Sa = T0 / 2 / T * Vdc;
        Sb = (T - T0 / 2) / T * Vdc;
        Sc = (T - T0 / 2 - T1) / T * Vdc;
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
        
        Sa = T0 / 2 / T * Vdc;
        Sc = (T - T0 / 2) / T * Vdc;
        Sb = (T - T0 / 2 - T2) / T * Vdc;
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
        
        C(1 : T0 / 4) = 0;
        C(T0 / 4 + 1 : (T0 / 4 + T1 / 2)) = 2 / 3 * Vdc;
        C(T0 / 4 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2) = 1 / 3 * Vdc;
        //C(T0 / 4 + T1 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2) = 0
        C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2) = 1 / 3 * Vdc;
        C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2) = 2 / 3 * Vdc;
        //C(T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + 1 : T0 / 4 + T1 / 2 + T2 / 2 + T0 / 2 + T2 / 2 + T1 / 2 + T0 / 4) = 0;
        
        Sb = T0 / 2 / T * Vdc;
        Sc = (T - T0 / 2) / T * Vdc;
        Sa = (T - T0 / 2 - T1) / T * Vdc;
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
        
        Sb = T0 / 2 / T * Vdc;
        Sa = (T - T0 / 2) / T * Vdc;
        Sc = (T - T0 / 2 - T2) / T * Vdc;
    end
        
    f1 = scf(1);
    plot(start_time : 1 / scale : start_time + T / scale, a + 2, 'r');
    plot(start_time : 1 / scale : start_time + T / scale, b, 'g');
    plot(start_time : 1 / scale : start_time + T / scale, c - 2, 'b');
    xtitle('svpwm curve', 'time', 'value');
    legend('phase a', 'phase b', 'phase c');
    
    f2 = scf(2);
    plot(start_time : 1 / scale : start_time + T / scale, A + 2 * Vdc, 'r');
    plot(start_time : 1 / scale : start_time + T / scale, B, 'g');
    plot(start_time : 1 / scale : start_time + T / scale, C - 2 * Vdc, 'b');
    xtitle('phase voltage', 'time', 'value');
    legend('phase a', 'phase b', 'phase c');
    //plot(start_time : 1 / scale : start_time + T / scale, N * n, '*c');
    
    f3 = scf(3);
    pSa = zeros(1 : 1 + start_time);
    pSa($) = Sa;
    plot(pSa, 'r*');
    pSb = zeros(1 : 1 + start_time);
    pSb($) = Sb;
    plot(pSb, 'g*');
    pSc = zeros(1 : 1 + start_time);
    pSc($) = Sc;
    plot(pSc, 'b*');
endfunction


