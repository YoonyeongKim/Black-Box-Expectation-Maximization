function [x] = RK4_thrust(Func, x, u1, u2, mass, delt)
k1  = Func(x, u1, u2, mass)*delt;
k2  = Func(x + k1*0.5, u1, u2, mass)*delt;
k3  = Func(x + k2*0.5, u1, u2, mass)*delt;
k4  = Func(x + k3, u1, u2, mass)*delt;
x = x + (k1 + 2*(k2 + k3) + k4)/6.0;
end

