global A;
global R;
global b;
tic
for theta=1:360
    B=A(:,theta)*pinv((A(:,theta))'*A(:,theta))*A(:,theta)';
    P_A_theta{theta}=B;
end
Y=trace((A*pinv(A'*A)*A')*R);
for i=1:360
    T=log(1+i);
    fun(i)=trace(P_A_theta{i}*R);
end   
j=1:1:360;
plot(1:360,fun,'-r','linewidth',2);
toc
