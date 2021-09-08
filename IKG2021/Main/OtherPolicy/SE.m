function [a_new, S_new, tau, Ybar] = SE(a_old, S, tau, Ybar, Ybar_update, T)

gamma = 1;

if a_old == S(end) % update Ybar for new loop
    Ybar = Ybar_update;
    tau = tau + 1;
    a_old = 0;
end
Ybar_max = max(Ybar(S));
S_new = S;

stop_flag = 0;
while 1
    Ybar_max_adj = Ybar_max - gamma * UtauT(tau,T);
    for i = S
        if i > a_old
            if Ybar(i) >= Ybar_max_adj
                a_new = i;
                stop_flag = 1;
                break;
            else
                S_new((S_new==i)) = [];
            end
        end
    end
    if stop_flag == 1
        break;
    else
        a_old = 0;
        tau = tau + 1;
        S = S_new;
    end
end

end


function U = UtauT (tau,T)
    a = T/tau;
    if a > exp(1)
        b = log(a);
    else
        b = 1;
    end
    U = 2 * sqrt(2*b/tau);
end