function [ real,imag ] = read_from_dat( data_name )
file_id=fopen(data_name,'rb');%打开这个iq文件
[row,count] = fread(file_id,inf,'short'); %
fclose(file_id);
I = row(1:2:count); 
Q = row(2:2:count);  
if length(I)>length(Q)
   I = I(1:length(I)-1); 
end

real = I';
imag = Q';


%C = I + 1i*Q;
% real = I/mean(abs(C));
% imag = Q/mean(abs(C));
% real = real';
% imag = imag';


end

