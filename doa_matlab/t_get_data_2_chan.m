function [X,R] = t_get_data_2_chan(Dir,points,fn)
%读取数据文件 双通道
   channels = 2;
   for i=1:channels
      [Dir,'\fft',num2str(i),'.dat'];
      [I{i},Q{i}]=read_from_dat([Dir,'\fft',num2str(i),'.dat']);  
   end
   %save('I','I');  

   thw_num = 8;
   

%=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x

%==================  校正 =================================================

  for k = 1:thw_num
       range = (k-1)*points+1 : k*points;
       for i=1:channels           
              I1_standard{i} = I{i}(range);
              Q1_standard{i} = Q{i}(range);
              IQ_standard{i} = I1_standard{i}+1i*Q1_standard{i};
       end            
       for i=1:channels
            standard = wise_dot_division(IQ_standard{i},IQ_standard{1});           
            mean_val(k,i) = mean(standard,2);
       end 
  end
  for i = 1:channels
       mean_vector_temp = [];
       for k = 1:thw_num
           mean_vector_temp = [mean_vector_temp, mean_val(k,i)*ones(1,points)];
       end
       mean_vector{i} = mean_vector_temp;
  end

  
  FRL1 = points*(thw_num);%
  FRL2 = points*thw_num;
  for i=1:channels
      I1{i} = I{i}(FRL1+1+(fn-1)*FRL1:FRL1+(fn-1)*FRL1+FRL2);
      Q1{i} = Q{i}(FRL1+1+(fn-1)*FRL1:FRL1+(fn-1)*FRL1+FRL2);
      IQ{i} = I1{i}+1i*Q1{i};
  end
  
%   figure(100)
%   plot(real(IQ{1}));
%   
%   %fn;%帧数
%   figure(101)
%   plot(imag(IQ{1}));
  
  X = [];
  for i=1:channels
        X = [X;IQ{i}./mean_vector{i}];
  end
  
  

%==================  测向 =================================================

%============================== 换算至频域 =================================

 
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
 showlevel = 20*log10(sqrt(sumall));
 M2=sumall;
 M4=sumall2;
 ratio=sqrt(2*M2^2-M4)/(M2-sqrt(2*M2^2-M4));
 SNR=10*log10(ratio);
 X1 = fft(XX1,2^nextpow2(points),2);
 X2 = fft(XX2,2^nextpow2(points),2);
 
%============================== 重构R矩阵 ==================================
for kk=1:points

    XXNew(1,kk)=1+0*1i;
    XXNew(2,kk)=X2(2,kk)/X1(2,kk);
    XXNew(3,kk)=X2(3,kk)/X1(3,kk);
    XXNew(4,kk)=X2(4,kk)/X1(4,kk);
    XXNew(5,kk)=X2(5,kk)/X1(5,kk);
    XXNew(6,kk)=X2(6,kk)/X1(6,kk);
    XXNew(7,kk)=X2(7,kk)/X1(7,kk);
    XXNew(8,kk)=X2(8,kk)/X1(8,kk);
    XXNew(9,kk)=X2(9,kk)/X1(9,kk);

    R1(:,:,kk) = XXNew(:,kk)*XXNew(:,kk)'*abs(X1(2,kk))^2;%

end
    R = sum(R1,3); 
end