%读取数据文件 九通道
selection{1}='500MHZ 103单点';
selection{2}='500MHZ 258.2单点';
selection{3}='500MHZ 318.8单点';
selection{4}='500MHZ 三个信号（10KHZ）';%2500次
selection{5}='500MHZ 三个信号（10KHZ）317,260相干';%800次
selection{6}='500MHZ 三个信号（30KHZ）';%1200次
selection{7}='500MHZ 两个信号（30KHZ）103,260';%700次
selection{8}='500MHZ 两个信号（10KHZ）317,260相干';%1000次
Dir='D:\桌面\暑假华日\MCMC-DOA\干涉仪MATLAB代码\九通道多信号数据\';
Dir=[Dir,selection{7}];
points=8192;
fn=2;
for i=1:9    
    [I{i},Q{i}]=read_from_dat([Dir,'\fft',num2str(i),'.dat']);
end

for i=1:9
    I1_standard{i} = I{i}(1:points);
    Q1_standard{i} = Q{i}(1:points);
    IQ_standard{i} = I1_standard{i} + 1i*Q1_standard{i};
end
for i=1:9
    standard = IQ_standard{i}./IQ_standard{1};
    mean_val(i) = mean(standard,2);
end
for i=1:9
    I1{i} = I{i}(fn*points+1:(fn+1)*points);
    Q1{i} = Q{i}(fn*points+1:(fn+1)*points);
    IQ{i} = I1{i}+1i*Q1{i};
end

X = [];
for i=1:9
 X = [X;IQ{i}/mean_val(i)];%9*8192
end
  
  

%==================  测向 =================================================

%============================== 换算至频域 =================================
 thw_num=1;
 for i=1:thw_num
     XX1(i+1,:) = X(1,(i-1)*points+1:i*points);
     XX2(i+1,:) = X(2,(i-1)*points+1:i*points);%
 end
 sumx1=zeros(1,thw_num);
 sumx2=zeros(1,thw_num);
 sumx3=zeros(1,thw_num);
 sumx4=zeros(1,thw_num);
 for i=1:thw_num
     for j=1:points
         sumx1(i)=sumx1(i)+abs(XX1(i+1,j))^2;
         sumx2(i)=sumx2(i)+abs(XX2(i+1,j))^2;
         sumx3(i)=sumx3(i)+abs(XX2(i+1,j))^4;
         sumx4(i)=sumx4(i)+abs(XX2(i+1,j))^4;
     end
 end
 sumxx1=sum(sumx1./points)/thw_num;
 sumxx2=sum(sumx2./points);
 
 sumxx3=sum(sumx3./points)/thw_num;
 sumxx4=sum(sumx4./points);
 
 sumall=(sumxx1+sumxx2)/(thw_num+1);  
 sumall2=(sumxx3+sumxx4)/(thw_num+1);
 showlevel = 20*log10(sqrt(sumall))
