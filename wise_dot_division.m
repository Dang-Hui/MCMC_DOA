function result = wise_dot_division(a,b)
result = [];
for i = 1:length(b)
    if b(i)~=0
        result = [result,a(i)/b(i)];
    end
end