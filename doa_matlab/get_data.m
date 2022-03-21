function [R,A]=get_data2(Dir_me)
%% 构建A矩阵，构建X矩阵
k = (0:8)';  
fn=2;% 帧数
points = 8192;   %校正
K=8192;
f = 500e6; % 500MHz
r = 0.58;
c = 2.998*10^8; % 光速
A_=@(reference)exp(1i*(2*pi*f*r/c*cos(2*pi*k/9-reference*pi/180)));%构建A
A=[];
for i=1:360
    A=[A,A_(i)];
end
save('A','A');
directory = Dir_me;
for i=1:9    
    [I{i},Q{i}]=read_from_dat([directory,'\fft',num2str(i),'.dat']);
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

X1 = [];
 for i=1:9
     X1 = [X1;IQ{i}/mean_val(i)];%9*8192
 end
R=X1*X1';
end