% TODO: make a modular function out of this

x0 = 0;
xf = 10;
v0 = 0;
vf = 0;
a0 = 0;
af = 0;
j0 = 0;
jf = 0;
ti = 0;
T = 50;


A = [0 0 0 0 0 0 0 1;
    T^7 T^6 T^5 T^4 T^3 T^2 T 1;
    0 0 0 0 0 0 1 0;
    7*T^6 6*T^5 5*T^4 4*T^3 3*T^2 2*T 1 0
    0 0 0 0 0 2 0 0;
    42*T^5 30*T^4 20*T^3 12*T^2 6*T 2 0 0;
    0 0 0 0 6 0 0 0;
    210*T^4 120*T^3 60*T^2 24*T 6 0 0 0];
    
    

B = [x0; xf; v0; vf; a0; af; j0; jf];

X = linsolve(A,B);

posTime = [X(1,1) X(2,1) X(3,1) X(4,1) X(5,1) X(6,1) X(7,1) X(8,1)];
velTime = [7*X(1,1) 6*X(2,1) 5*X(3,1) 4*X(4,1) 3*X(5,1) 2*X(6,1) X(7,1)];
accelTime = [42*X(1,1) 30*X(2,1) 20*X(3,1) 12*X(4,1) 6*X(5,1) 2*X(6,1)];
jerkTime = [210*X(1,1) 120*X(2,1) 60*X(3,1) 24*X(4,1) 6*X(5,1)];

subplot(4,1,1);
x = linspace(ti,T);
y1 = posTime;
plot(x, polyval(y1,x))
title("position vs time")
xlabel("t")
ylabel("x(t)")


subplot(4,1,2);
y2 = velTime;
plot(x, polyval(y2,x))
title("velocity vs time")
xlabel("t")
ylabel("x'(t)")

subplot(4,1,3);
y3 = accelTime;
plot(x, polyval(y3,x))
title("acceleration vs time")
xlabel("t")
ylabel("x''(t)")

subplot(4,1,4);
y4 = jerkTime;
plot(x, polyval(y4,x))
title("jerk vs time")
xlabel("t")
ylabel("x'''(t)")
