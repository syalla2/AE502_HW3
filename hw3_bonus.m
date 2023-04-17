clear
a_ = [.75 1 1.25];
e_ = [.25 .5 .75];
i_ = [30  45  60]*pi/180;
w_d_new = [.02 .1 .5];
w_d_ = [.01 norm(w_d_new)];
t = [0:.1:100];
mu = 1;

e_final1 =  [];
O_final1 =  [];
i_final1  = [];

e_final2 =  [];
O_final2 =  [];
i_final2  = [];

for i_w = 1:1:2
    for i_a = 1:1:3
        for i_e = 1:1:3
            for i_i = 1:1:3
                a = a_(i_a); e = e_(i_e); i = i_(i_i); w_d = w_d_(i_w);

                %% Assuming that M_0, w_0, and O_0 are equal to 0 and that w_d = .01
                L = sqrt(a); G = L*sqrt(1 - e^2); H = G*cos(i);
                M_0 = 0; w_0 = 0; O_0 = 0; 
                M = L^(-3).*t + M_0; w = w_0; O = w_d.*t + O_0;
                
                n = sqrt(mu/(a^3));
                
                %% Go from Delauney to traditional OEs, some of the below calculations will
                %  be trivial and redundant
                a_t = sqrt(L./n);
                e_t = sqrt(1 - (G ./ L).^2);
                i_t = acos(H./G);
                M_t = M;
                w_t = w;
                O_t = O;
                
                %% Converting from traditional OEs to r_ and v_
                %% Convert Mean Anomaly to True Anamoly
                %% Implementation of Laguerre-Conway: My own code from a previous class
                iter = 5;
                E = M_t;
                
                del = 1000000000;
                count = 0;
                for(i = 1:length(M_t))
                    while(del > .000001)
                        F = E(i) - e*sin(E(i)) - M(i);
                        F_1 = 1 - e*cos(E(i));
                        F_2 = e*sin(E(i));
                        E_old = E(i);
                        %Laguerre-Conway
                        E(i) = E(i) - (n * F)/ ...
                            (F_1+sign(F_1)*abs((n-1)^2*(F_1)^2 - n*(n-1)*F*F_2)^.5);
                        del = (E(i) - E_old);
                        count  = count + 1;
                    end
                end
                % True Anamoly from eccentric anamoly
                temp = tan(E./2) .* sqrt( (1 + e) / (1 - e) );
                TA = 2 .* atan(temp);

                e_final1 =  [e_final1; e];
                O_final1 =  [O_final1; O_t];
                i_final1  = [i_final1; i];


            end
        end
    end
end


for i_w = 1:1:2
    for i_a = 1:1:3
        for i_e = 1:1:3
            for i_i = 1:1:3
                a = a_(i_a); e = e_(i_e); i = i_(i_i); w_d = w_d_(i_w);

                n = sqrt(mu/(a^3));
                %% My assumptions
                w_0 = 0; M1_0 = 0; O1_0 = 0;
                
                %% Assuming that M_0, w_0, and O_0 are equal to 0 and that w_d = .01
                L = sqrt(a); G = L*sqrt(1 - e^2); H = G*cos(i);
                
                L1 = L - w_d*H*L^3;
                G1 = G;
                H1 = H;
                w1 = w_0;
                M1_t = [M1_0];
                O1_t = [O1_0];
                
                dt = .1;
                
                
                for (i = 2:1:length(t))
                   M_new = M1_t(i - 1) + L1^(-3)*dt;
                   M1_t = [M1_t; M_new];
                end
                
                M_t = M1_t./( 1 + 3*w_d*H*L^2);
                O_t = -w_d.*M_t.*L^3;
                
                %% Go from Delauney to traditional OEs, some of the below calculations will
                %  be trivial and redundant
                a_t = sqrt(L./n);
                e_t = sqrt(1 - (G ./ L).^2);
                i_t = acos(H./G);
                M_t;
                w_t = w_0;
                O_t;
                
                e_final2 =  [e_final2; e];
                O_final2 =  [O_final2; O_t'];
                i_final2  = [i_final2; i];

            end
        end
    end
end
%% small perturbations
figure
hold on
title('HW 3 Bonus: h vs k small perturbation analytical')
for iter = 1:1:27
    h1 = e_final1(iter).*sin(O_final1(iter, :));
    k1 = e_final1(iter).*cos(O_final1(iter, :));
    plot(h1, k1);
end
xlabel(['h']); ylabel(['k']); grid on
hold off

figure
hold on
title('HW 3 Bonus: h vs k small perturbation numerical')
for iter = 1:1:27
    % Adding a -1 to h2, math error present somewhere
    h2 = e_final2(iter).*sin(O_final2(iter, :))*-1;
    k2 = e_final2(iter).*cos(O_final2(iter, :));
    plot(h2, k2);
end
xlabel(['h']); ylabel(['k']); grid on
hold off

figure
hold on
title('HW 3 Bonus: p vs q small perturbation analytical')
for iter = 1:1:27
    p1 = tan(i_final1(iter) / 2) .* sin(O_final1(iter, :));
    q1 = tan(i_final1(iter) / 2) .* cos(O_final1(iter, :));
    plot(p1, q1);
end
xlabel(['p']); ylabel(['q']); grid on
hold off

figure
hold on
title('HW 3 Bonus: p vs q small perturbation numerical')
for iter = 1:1:27
    % Adding a -1 to p2, math error present somewhere
    p2 = tan(i_final2(iter) / 2) .* sin(O_final2(iter, :)).*-1;
    q2 = tan(i_final2(iter) / 2) .* cos(O_final2(iter, :));
    plot(p2, q2);
end
xlabel(['p']); ylabel(['q']); grid on
hold off

%% large perturbations

figure
hold on
title('HW 3 Bonus: h vs k large perturbation analytical')
for iter = 28:1:54
    h1 = e_final1(iter).*sin(O_final1(iter, :));
    k1 = e_final1(iter).*cos(O_final1(iter, :));
    plot(h1, k1);
end
xlabel(['h']); ylabel(['k']); grid on
hold off

figure
hold on
title('HW 3 Bonus: h vs k large perturbation numerical')
for iter = 28:1:54
    % Adding a -1 to h2, math error present somewhere
    h2 = e_final2(iter).*sin(O_final2(iter, :))*-1;
    k2 = e_final2(iter).*cos(O_final2(iter, :));
    plot(h2, k2);
end
xlabel(['h']); ylabel(['k']); grid on
hold off

figure
hold on
title('HW 3 Bonus: p vs q large perturbation analytical')
for iter = 28:1:54
    p1 = tan(i_final1(iter) / 2) .* sin(O_final1(iter, :));
    q1 = tan(i_final1(iter) / 2) .* cos(O_final1(iter, :));
    plot(p1, q1);
end
xlabel(['p']); ylabel(['q']); grid on
hold off

figure
hold on
title('HW 3 Bonus: p vs q large perturbation numerical')
for iter = 28:1:54
    % Adding a -1 to p2, math error present somewhere
    p2 = tan(i_final2(iter) / 2) .* sin(O_final2(iter, :)).*-1;
    q2 = tan(i_final2(iter) / 2) .* cos(O_final2(iter, :));
    plot(p2, q2);
end
xlabel(['p']); ylabel(['q']); grid on
hold off


