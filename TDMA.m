 function X = TDMA(a,b,c,d)
n=length(d);
c(1)=c(1)/b(1);
d(1)=d(1)/b(1);
% for i=2:n-1
%     c(i)=c(i)/(b(i)-a(i)*c(i-1));
% end

for i=2:n
    c(i)=c(i)/(b(i)-a(i)*c(i-1));
    d(i)=(d(i)-a(i)*d(i-1))/(b(i)-a(i)*c(i-1));
end
X(n)=d(n);
for i=n-1:-1:1
    X(i)=d(i)-X(i+1)*c(i);
end