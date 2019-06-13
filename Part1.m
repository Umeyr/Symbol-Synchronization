

step = 0.01;
T = 1;
alfa = 1;
k = 50 ;
t = 0:step:T-step;

%pulse
p = alfa*exp(-k*((t-0.5*T).^2));

%information bits
M = 100;
a = rand(1, M);
a = 1 - 2*round(a);


s = zeros(1,M*size(t,2));
%generation of s(t)
for i = 1:M

s((i-1)*size(t,2)+1:(i-1)*size(t,2)+size(t,2)) =   a(i) * p;

end

%plot for first ten symbols
figure(1);
xaxis = [0:.01:10-.01];
plot(xaxis,s(1,1:1000));
xlabel('time');
ylabel('symbols');
title('Plot for first 10 symbols');

MSE = zeros(1,100);
N0 = zeros(1,100);

for w = 1:100   % this for loop is used for MSE vs N0 and M
    
MSE(w) = 0;
for q = 1:1000

%generation of r(t) with addition of noise

 tao =   0.25*rand*T ; %delay
 
 tao = round(tao,2);
 

r = zeros (1, size(s,2) + tao*(T/step));
n = randn (1, size(s,2) + tao*(T/step));



r(1, tao*(T/step)+1: size(s,2) + tao*(T/step)) = s(1,:);

Es = sum(r.^2) / (size(s,2) + tao*(T/step));

%SNR = rand*100;

%N0 = Es/SNR;

N0(w) = 0.1*w ;

n = sqrt(N0(w)/2).*n;

r = r + n;

%since p is symmetric, impulse response of match filter is the same as p

y = conv(p,r);

%early late gate synchronization

e = 0.05;
delta = 0.01;
x =  -.25 + 0.5*rand*T ; %x is between -0.25 T and 0.25 T

t1 = round( x + T) ;
t2 = t1 + e;


for i = 1:M
    
    
if abs(y(round(t1*size(t,2),1))) > abs(y(round(t2*size(t,2),1)))

    t1 = t1 - delta + T;
    t2 = t1 + e;
elseif abs(y(round(t1*size(t,2),1))) < abs(y(round(t2*size(t,2),1)))
    
    t1 = t1 + delta + T;
    t2 = t1 + e;
else
     t1 = t1 + T;
     t2 = t1 + e;
    break;
end

end

tdelay = abs(((t1 + t2)/2 - (i+1)*T - 0.01));

MSE(w) = MSE(w) + (tdelay-tao)^2;

end

MSE(w) = MSE(w) / 1000;

end

 figure(2);
 plot(N0,MSE);
 xlabel('N0');
 ylabel('MSE');
 title('MSE vs N0');

