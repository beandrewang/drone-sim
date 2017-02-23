// implement a motor

exec('../inverter/inverter.sci');

function [Ia, Ib, Ic, M, theta] = pmsm_motor(A, B, C, R)
    // this is a simulator for a pmsm
    // A, B, C, the 3 phase voltages 
    // R, the sample resistance.
    // Ia, Ib, Ic, the 3 phase current
    // M, the resultance voltage vector magnitude 
    // theta, the resultance voltage vector phase
    
    [Ua, Ub] = forward_clark_transform(A, B, C);
    Magnitude = sqrt(Ua ^ 2 + Ub ^ 2);
    theta = atan(Ub, Ua);
    
    Ia = Ua / R;
    Ib = Ub / R;
    Ic = Uc / R;
endfunction
