function[err] = err_norm(yo,y1,y1_hat,Atol,Rtol)
n = length(yo);
sc = zeros(n,1);

% implement sci 
for i=1:n
    sc(i) = Atol + max(abs(yo(i)),abs(y1(i)))*Rtol;
end

% compute err
err = sqrt (1/n * sum(((y1 - y1_hat)./sc).^2));

end
