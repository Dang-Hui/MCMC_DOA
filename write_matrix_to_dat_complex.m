function write_matrix_to_dat_complex(data_name,vari)
x=[];
for i=1:length(vari(:,1))
    x=[x,vari(i,:)]
end
vari=x;
file_id=fopen([data_name,'.dat'],'wb');
Long=length(vari);
A(1:2:2*Long)=real(vari)';
A(2:2:2*Long)=imag(vari)';
fwrite(file_id,A,'double');
fclose(file_id);
end

