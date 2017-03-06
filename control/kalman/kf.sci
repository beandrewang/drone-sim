// this is a implement of a kalman filter

function kf = kalman_filter_init(A, B, C, X0, P0, u, Q, R)
    // A, the system matrix
    // B, the control matrix
    // C, the transfer matrix
    // X0, the initial state of the system
    // P0, the coinvariance matrix of the initial state, default 0.
    // u, the control vector of the system
    // Q, the coinvariant of the predict, there might be some disturbance or errors 
    // R, the coinvariant of the measured value
    // H, the scale matrix between the state and the measured value
    // kf, the struct of the kalman filter
    kf = struct('A', A, 'B', B, 'C', C, 'X', X0, 'u', u, 'P', P0, 'Q', Q, 'R', R);
endfunction

function kf = kalman_filter_predict(kf)
    // kf, the struct of the kalman filter
    kf.X = kf.A * kf.X + kf.B * kf.u;
    kf.P = A * kf.P * A' + kf.Q;
endfunction

function kf = kalman_filter_update(kf, Z)
    // kf, the struct of the kalman filter
    // K, the kalman gain
    // Z, the measured vector
    K = kf.P * kf.C' * inv( kf.C * kf.P * kf.C' + kf.R);
    kf.X = kf.X + K * (Z - kf.C * kf.X);
    kf.P = kf.P - K * kf.C * kf.P;
endfunction


