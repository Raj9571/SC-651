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
A=V*W';
[U,S,Z]=svd(A);%there is no direct way to compute polar decompostion in matlab
% so we will use SVD to compute it
U_pol=Z*U';
P_pol=U*S*U';%A=U_pol*P_pol

if det(U_pol)>0
    X=eye(3,3);
else
    X=[1 0 0; 
        0 1 0;
        0 0 -1];
end
 M=U*X*U'*U_pol;
 W_pred=M*V;

loss=0;
for i=1:6
    loss=loss+norm(W(:,i)-M*V(:,i))^2;
end
M
loss_whaba=loss