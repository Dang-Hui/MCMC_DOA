a=max(svd(R));
j=0;
while(abs(a/10)>1)
    j=j+1;
    a=a/10;
end
b=10^(j-1);
for theta=1:360
    P_A_theta=A(:,theta)*pinv(A(:,theta)'*A(:,theta))*A(:,theta)';
    p_theta(theta)=exp(1*trace(P_A_theta*R)/b);
end
max_a=max(p_theta);
min_a=min(p_theta);
mean_a=mean(p_theta);
plot(1:360,(p_theta-mean_a)/(max_a-min_a));
