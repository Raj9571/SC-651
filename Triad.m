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
% All vectors are unit vectors now
V=V_norm;
W=W_norm;

v1=[0.63, 0.45, 0.87];
v2=[0.03, 0.40, 0.86];

r1=v1/norm(v1);
r2=cross(v1,v2)/(norm(cross(v1,v2)));

r3=cross(r1,r2);

M_ref_T=[r1;r2;r3];


w1=[0.17, 0.03, 0.50];
w2=[0.26, 0.37, 0.15];

s1=w1/norm(w1);
s2=cross(w1,w2)/(norm(cross(w1,w2)));
s3=cross(s1,s2);

M_obs_T=[s1;s2;s3];
M_obs=M_obs_T.';
A=M_obs*M_ref_T;


loss=0;
for i=1:6
loss=loss+norm(W(:,i)-A*V(:,i))^2;
end

A
MSE_triad=loss