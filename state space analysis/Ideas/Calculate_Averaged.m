A = [1 2; 3 4];
B = [5 6]';
C = [0 0; 0 0];
D = 0;
u = 5;
x0 = [7 9]';
T = 1;%1e-3;
npts = 10000;
t = linspace(0,T,npts);
[y, t, x] = lsim(ss(A,B,C,D), u*ones(1,npts), t, x0);

avgx = mean(x,1)'

avgx2 = 1/T*((inv(A)*expm(A*T)-inv(A))*x0 + inv(A)*(inv(A)*(expm(A*T)-eye(2))-T*eye(2))*B*u)

avgx3 = inv(A)/T*((expm(A*T)-eye(2))*x0 + (inv(A)*(expm(A*T)-eye(2))-T*eye(2))*B*u)

avgx4 = (1/T)*A\(    (expm(A*T)-eye(2))*x0 +    A\(expm(A*T)-eye(2))*B*u -     T*eye(2)*B*u   )

avgx5 = (1/T)*A\(  (expm(A*T)-eye(2))*(x0 + A\B*u) -     T*B*u   )

T*A^2*avgx5 - (  (expm(A*T)-eye(2))*(A*x0 + B*u) -     A*T*B*u   )

[X,flag] = lsqr(T*A^2, (  (expm(A*T)-eye(2))*(A*x0 + B*u) -     A*T*B*u   ))

T*A^2*X - (  (expm(A*T)-eye(2))*(A*x0 + B*u) -     A*T*B*u   )



diff = avgx4 - avgx%T*A*avgx - ((expm(A*T)-eye(2))*x0 + (inv(A)*(expm(A*T)-eye(2))-T*eye(2))*B*u)



% avgx4 = A*x +B*u



