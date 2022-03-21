%% 读取数据文件
selection{1}='500MHZ 103单点';
selection{2}='500MHZ 258.2单点';
selection{3}='500MHZ 318.8单点';
selection{4}='500MHZ 三个信号（10KHZ）';%2500次
selection{5}='500MHZ 三个信号（10KHZ）317,260相干';%800次
selection{6}='500MHZ 三个信号（30KHZ）';%1200次
selection{7}='500MHZ 两个信号（30KHZ）103,260';%700次
selection{8}='500MHZ 两个信号（10KHZ）317,260相干';%1000次
Dir_me='D:\桌面\暑假华日\MCMC-DOA\干涉仪MATLAB代码\九通道多信号数据\';
Dir_me=[Dir_me,selection{7}];
global R;
global A;
% global TOT mchain;
[R,A]=get_data(Dir_me);
% [X,R] = t_get_data_2_chan(Dir_me,8192,2);
%% 确定模型参数
% a=max(max(R));
a=R(1,1);
% a=max(svd(R));
j=0;
while(abs(a/10)>=1)
    j=j+1;
    a=a/10;
end
global b;
b=10^j;
tic
direction=MCMC(2,1000);%设置信号源数和迭代次数
toc
%% MCMC算法
function direction= MCMC(sig_num,iteration)
start_value=ceil(rand(1,sig_num)*360);
flag=0;
burn_in=iteration-100;
while 1 && sig_num~=1
    mcmc_chain=mh(start_value,iteration,sig_num);
    for i=1:sig_num
        table=tabulate(mcmc_chain(burn_in:iteration,i));
    %     figure;
    %     plot(mcmc_chain(:,i),'-k','linewidth',1);
    %     hold on;
        [~,idx]=max(table(:,2));
        theta=table(idx);
        direction(i)=theta;
    end
    x=nchoosek(direction,2);
    for kk=1:size(x,1)
        if abs(x(kk,1)-x(kk,2))<10
            break;
        end
        if kk==size(x,1)
            flag=1;
        end
    end
    if flag==1
        break;
    end
end
if sig_num==1
    mcmc_chain=mh(start_value,iteration,sig_num);
    for i=1:sig_num
        table=tabulate(mcmc_chain(burn_in:iteration,i));
    %     figure;
    %     plot(mcmc_chain(:,i),'-k','linewidth',1);
    %     hold on;
        [~,idx]=max(table(:,2));
        theta=table(idx);
        direction(i)=theta;
    end
end
save('mcmc_chain','mcmc_chain');
% 
% % burn_in=400;
% direction=[];
% for i=1:sig_num
%     table=tabulate(mcmc_chain(burn_in:iteration,i));
%     [~,idx]=max(table(:,2));
%     theta=table(idx);
%     direction=[direction,theta];
%     save('direction','direction');
% end
% global TOT mchain;
% mcmc_chain=mchain;
mcmc_chain(end,:)
for c=1:sig_num
    h(c)=plot(mcmc_chain(:,c),'-k','linewidth',1);
    hold on;
end
F=plot([1,1500],[317,317],'--k','linewidth',1);
 hold on;
plot([1,1500],[103,103],'--k','linewidth',1);
hold on;
plot([1,1500],[260,260],'--k','linewidth',1);
hold on;
legend([h(1),F],'估计值','真实值');
xlabel('样本点');
ylabel('方位角\(°)');
% title('MCMC三个信号源方位角估计');
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
global A;
global R;
% save('R','R');
% A_theta=[];
% for k=1:sig_num
%     A_theta=[A_theta,A(:,theta(k))];
% end
% save('A_theta','A_theta');
% P_A_theta=A_theta*pinv((A_theta)'*A_theta)*A_theta';   
global b;
% p_theta=exp(2*trace(P_A_theta*R)/b);
% save('p_theta','p_theta');
for k=1:sig_num
%     A_theta=[A_theta,A(:,theta(k))];
    P_A_theta{k}=A(:,theta(k)) * pinv(A(:,theta(k))'*A(:,theta(k))) * A(:,theta(k))';
    p_theta(k)=exp(5*trace(P_A_theta{k}*R)/b);
end
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
%     lanmuda=target_fun(proposal,m,sig_num)/target_fun(chain(m,:),m,sig_num);%接受比率
    p_theta1=target_fun(proposal,m,sig_num);
    p_theta2=target_fun(chain(m,:),m,sig_num);
    lan=1;
    for jj=1:sig_num
        lanmuda(jj)=p_theta1(jj)/p_theta2(jj);
        lan = lan*lanmuda(jj);
    end
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
