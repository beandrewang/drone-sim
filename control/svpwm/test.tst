// unite test for svpwm_generator
funcprot(0);
exec('svpwm_generator.sci');

theta = -%pi / 3;
Vd = 0;
Vq = 8;
T = 0.2;
Vdc = 12;
t = 0;

//for i = 0 : 1
    for theta = 0 : 0.1 : 2 * %pi
        [_alpha, _beta] = reverse_park_transform(Vd, Vq, theta);
        [a, b, c, t0, t1, t2, n, clock_freq] = svpwm_period_generator(_alpha, _beta, T, Vdc);
        [A, B, C] = plot_svpwm(t0, t1, t2, n, Vdc, t);
        mprintf('%d', n);
        t = t + T;
    end
//end




