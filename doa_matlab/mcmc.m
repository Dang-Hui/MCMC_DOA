function [direction,amplitudes]= mcmc(R,A,sig_num)
%确定模型参数
a=max(svd(R));
j=0;
while(abs(a/10)>1)
    j=j+1;
    a=a/10;
end
b=10^(j-1);
iteration= 200*sig_num^3;
start_value=ceil(rand(1,sig_num)*360);
mcmc_chain=mh(start_value,iteration,sig_num);
burn_in=iteration-100;
for i=1:sig_num
    table=tabulate(mcmc_chain(burn_in:iteration,i));
    [~,idx]=max(table(:,2));
    theta=table(idx);
    direction(i)=theta;
end
for ii=1:sig_num
    a=A(:,direction(ii));
    Q=pinv(a'*a)*a';
    S=real(Q*R*Q');
    amplitudes(ii)=S;
end
%% 提议分布
%自适应随机游走采样方法
function proposal = proposalfunction(theta,i,iteration,sig_num)
sig_max=180;
sig=sig_max*exp(2*log(sig_max)*(((iteration-i)/iteration)-1));
SIG=ones(1,sig_num)*sig;
theta2=mod(abs(round(normrnd(theta,SIG))),361);
while ~all(theta2)
    theta2=mod(abs(round(normrnd(theta,SIG))),361);
end
proposal=theta2;
end
%独立马尔科夫链采样方法
function proposal=proposalfunction2(~,~,~,sig_num)
theta2=abs(round(unifrnd(1,360,1,sig_num)));
proposal=theta2;
end
%% 目标函数
function p_theta=target_fun(theta,~,sig_num)
A_theta=[];
for k=1:sig_num
    A_theta=[A_theta,A(:,theta(k))];
end
save('A_theta','A_theta');
P_A_theta=A_theta*pinv((A_theta)'*A_theta)*A_theta';   
p_theta=exp(5*trace(P_A_theta*R)/b);
save('p_theta','p_theta');
end

%% M-H
function mcmc_chain=mh(startvalue,iterations,sig_num)
chain=zeros(iterations,sig_num);
chain(1,:)=startvalue;
for m=1:iterations
    u2=rand();
    if u2<0.5
        proposal = proposalfunction2(chain(m,:),m,iterations,sig_num);
    else
        proposal = proposalfunction(chain(m,:),m,iterations,sig_num);
    end
    lanmuda=target_fun(proposal,m,sig_num)/target_fun(chain(m,:),m,sig_num);%接受比率
    alpha=min(1,lanmuda);%接受概率
    u=rand();
    if(u<=alpha) 
       chain(m+1,:)=proposal;
    else
       chain(m+1,:)=chain(m,:);
    end
end
mcmc_chain=chain();
end

end

