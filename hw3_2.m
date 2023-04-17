clear
%% Givens from problem
a = 1; e = .5; i = 45*pi/180; t = [0:.001:100]; mu = 1;

%% Assuming that M_0, w_0, and O_0 are equal to 0 and that w_d = .01
L = sqrt(a); G = L*sqrt(1 - e^2); H = G*cos(i);
M_0 = 0; w_0 = 0; O_0 = 0; w_d = .01;
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
mag_r = a*(1 - e.*cos(E));
h = sqrt(mu*a*(1 - e^2));
r_X = mag_r.*(cos(O_t).*cos(w_t + TA) - sin(O_t).*sin(w_t + TA)*cos(i));
r_Y = mag_r.*(sin(O_t).*cos(w_t + TA) + cos(O_t).*sin(w_t + TA)*cos(i));
r_Z = mag_r.*(sin(i) .* sin(w_t + TA));

figure
hold on
plot3(r_X, r_Y, r_Z)
xlabel('X [Canonical Distance Unit]'); ylabel('Y [Canonical Distance Unit]');
zlabel('Z [Canonical Distance Unit]'); grid on; title('Homework 3: Problem 2 Orbit in Cartesian Space')
hold off
