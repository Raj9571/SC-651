V_t=[0.63, 0.45, 0.87; 
    0.03, 0.40, 0.86;
    0.86, 0.17, 0.16;
    0.29, 0.79, 0.04; 
    0.60, 0.38, 0.15; 
    0.99, 0.43, 0.76];

W_t=[0.17, 0.03, 0.50;
    0.26, 0.37, 0.15;
    0.78, 0.65, 0.67; 
    0.83, 0.96, 0.82;
    0.69, 0.72, 0.04;
    0.64, 0.77, 0.34];
V=V_t';
W=W_t';

a=1/6;
V_norm=zeros(3,6);
for i=1:6
    V_norm(:,i)=V(:,i)/norm(V(:,i));
end
W_norm=zeros(3,6);
for i=1:6
    W_norm(:,i)=W(:,i)/norm(W(:,i));
end
V=V_norm;
W=W_norm;
B=a*W*V';
sigma=trace(B);
S=B+B';
Z=zeros(3,1);
for i=1:6
Z=Z+a*cross(W(:,i),V(:,i));
end

k=trace(adjoint(S));
delta=det(S);
a=sigma^2-k;
d=Z'*S*S*Z;
c=delta+Z'*S*Z;
b=sigma^2+Z'*Z;
syms x;
f=x^4-(a+b)*x^2-c*x+(a*b+c*sigma-d);

g=diff(f);
epsilon=10^(-5);
guess=1;
lambda=New_Raph(f,g,epsilon,guess);

 omega = lambda;
 alpha = omega^2 - sigma^2 + k;

    beta = omega - sigma;
    gamma = (omega + sigma) * alpha - delta;

    % (Eqn.68) Compute X vector 
    X = (alpha * eye(3) + beta * S + S^2) * Z;

    % Compute the optimal quaternion (Eqn.69)
    q = [X; gamma] ./ sqrt(gamma^2 + norm(X)^2);

    % Convert the optimal quaternion to DCM
    % standard method
    q0 = q(4);
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
 A=zeros(3,3);
 A(1,1) = q0^2 + q1^2 - q2^2 - q3^2;
    A(1,2) = 2*(q1*q2 + q0*q3);
    A(1,3) = 2*(q1*q3 - q0*q2);
    A(2,1) = 2*(q1*q2 - q0*q3);
    A(2,2) = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    A(2,3) = 2*(q2*q3 + q0*q1);
    A(3,1) = 2*(q1*q3 + q0*q2);
    A(3,2) = 2*(q2*q3 - q0*q1);
    A(3,3) = q0^2 - q1^2 - q2^2 + q3^2;

loss=0;
for i=1:6
loss=loss+norm(W(:,i)-A*V(:,i))^2;
end
A
MSE_quest=loss

function y=New_Raph(f,g,epsilon,guess)
syms x;
for i=1:100
     f0=vpa(subs(f,x,guess)); %Calculating the value of function at x0
     f0_der=vpa(subs(g,x,guess)); %Calculating the value of function derivative at x0
  y=guess-f0/f0_der; % The Formula
err=abs(y-guess);
if err<epsilon %checking the amount of error at each iteration
break
end
guess=y;
end
end

