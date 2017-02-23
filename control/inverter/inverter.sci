// implementation for a inverter

  
function [A, B, C] = l2p3_inverter(a, b, c, t0, t1, t2, n, Vdc)
    // this is a two level three phase VIS inverter
    // a, the 1st upper MOSFET control
    // b, the 2nd upper MOSFET control
    // c, the 3rd upper MOSFET control
    // t0, the 0 component time
    // t1, the 1 component time
    // t2, the 2 component time 
    // n, current sector
    // Vdc, the DC voltage source for the inverter
    // A, B, C, the 3 phase voltages
    
    clock_freq = (t0 + t1 + t2) / 1000;
    scale = 1 / clock_freq;
    T = ceil((t0 + t1 + t2) * scale);
    
    A = zeros(1, T + 1);
    B = zeros(1, T + 1);
    C = zeros(1, T + 1);
    
    T0 = t0 * scale;
    T1 = t1 * scale;
    T2 = t2 * scale;
    
    if 1 == n then
        //sector 1
        SA = (T1 / T * 2 / 3 + T2 / T * 1 / 3) * Vdc;
        SB = (-T1 / T * 1 / 3 + T2 / T * 1 / 3) * Vdc;
        SC = (-T1 / T * 1 / 3 - T2 / T * 2 / 3) * Vdc;
    elseif 2 == n then
        SA = (-T2 / T * 1 / 3 + T1 / T * 1 / 3) * Vdc;
        SB = (T2 / T * 2 / 3 + T1 / T * 1 / 3) * Vdc;
        SC = (-T2 / T * 1 / 3 - T1 / T * 2 / 3) * Vdc;
    elseif 3 == n then
        SA = (-T1 / T * 1 / 3 - T2 / T * 2 / 3) * Vdc;
        SB = (T1 / T * 2 / 3 + T2 / T * 1 / 3) * Vdc;
        SC = (-T1 / T * 1 / 3 + T2 / T * 1 / 3) * Vdc;
    elseif 4 == n then
        SA = (-T2 / T * 1 / 3 - T1 / T * 2 / 3) * Vdc;
        SB = (-T2 / T * 1 / 3 + T1 / T * 1 / 3) * Vdc;
        SC = (T2 / T * 2 / 3 + T1 / T * 1 / 3) * Vdc;
    elseif 5 == n then
        SA = (-T1 / T * 1 / 3 + T2 / T * 1 / 3) * Vdc;
        SB = (-T1 / T * 1 / 3 - T2 / T * 2 / 3) * Vdc;
        SC = (T1 / T * 2 / 3 + T2 / T * 1 / 3) * Vdc;
    elseif 6 == n then
        SA = (T2 / T * 2 / 3 + T1 / T * 1 / 3) * Vdc;
        SB = (-T2 / T * 1 / 3 - T1 / T * 2 / 3) * Vdc;
        SC = (-T2 / T * 1 / 3 + T1 / T * 1 / 3) * Vdc;
    end
    
endfunction
