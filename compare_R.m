I=20;          %ratio between Th and Tg

N=16;           %number of elements
M=8;            %number of antanas
Pt=10;        %transmit power
period=10000;
fade_var = 0.1; % variance of  Rayleigh channel


MSE=zeros(I,period);
MSE2=zeros(I,period);
MSE3=zeros(I,period);


sigma_W=sqrt(0.1);

sigma_V=0.05;
q=0.001;
h=zeros(N,N,I);
Q=sigma_V^2*eye(N);

phi = sqrtm(Pt) * dftmtx(16)/4;      %get the optimal phi based on the lemma 1
Y=zeros(M,N,I);

for i=1:N
    
    for j=1:M
        G(j,i) = sqrt(fade_var/2) * (normrnd(0,1,1,1) + 1i * normrnd(0,1,1,1)); 
        %G(j,i) =  (1 + 1i)*1/sqrt(2) ;        %Generate constant channel G
    end
end

for i=1:N
h1(i,1) = sqrt(0.5)/sqrt(2)* (normrnd(0,1,1,1) + 1i * normrnd(0,1,1,1));  %Generate random channel h

end


h(:,1,1)=h1;

for j=2:N
    
    h1=h1;
    h(:,j,1)=h1;
end

for b=2:I
    for j=1:N
    n=(1/sqrt(2))*normrnd(0,sigma_V,N,1)+ (1/sqrt(2))*1i * normrnd(0,sigma_V,N,1); % Inside noise
    
    h1=h1+h1.*n;
    
    h(:,j,b)=h1;
    
     W=(1/sqrt(2))*normrnd(0,sigma_W,M,1)+ (1/sqrt(2))*1i * normrnd(0,sigma_W,M,1); % Observation noise
     C(:,:,j)=G*diag(phi(j,:));
     Y1=C(:,:,j)*h(:,j,b)+W;
     Y(:,j,b)=Y1;
    
    end
end


for e=1:6
if e==1
    R=1;
elseif e==2
    R=0.75;
elseif e==3
    R=0.5;
elseif e==4
    R=0.25;
elseif e==5
    R=0.12;
elseif e==6
    R=0.1;

end


    
for x=1:period           % 5000 Monti Karlo

%Y=zeros(M,1);   %received signal

%C=zeros(M,N);



h_est=zeros(N,1);      %channel estimation

P_pri=zeros(N,N);      %prior covariance
P_post=zeros(N,N);      %posterior covariance
K=zeros(N,M);           %Kalman filter




% generate Rayleigh fading channel



 

%kalman filter algorithm
P_post=0.1*eye(16);
average_value = mean(h(:,1,1), 'all')/2;
h_est = zeros(16,1) * average_value;


MSE(1,x)=(norm(h(:,1,1)-h_est))^2;

for i=2:I
   
    for o=1:16
       
        Q=q*eye(N)*diag((h_est.*conj(h_est)));
         
        P_pri=P_post+Q;
        K=P_pri*(C(:,:,o))'*inv((C(:,:,o))*P_pri*(C(:,:,o))'+R*eye(M));
       
        h_est=h_est+K*(Y(:,o,i)-(C(:,:,o))*h_est);

        P_post=(eye(16)-K*(C(:,:,o)))*P_pri;
    
    end
    
    MSE(i,x)=(norm(h(:,16,i)-h_est,2))^2;
    



end
end
if e==1
    MSE1=mean(MSE,2);
elseif e==2
    MSE2=mean(MSE,2);
elseif e==3
    MSE3=mean(MSE,2);
elseif e==4
    MSE4=mean(MSE,2);
elseif e==5
    MSE5=mean(MSE,2);
elseif e==6
    MSE6=mean(MSE,2);

end
end
q1=1:1:I;
plot(q1,MSE1);
hold on
plot(q1,MSE2);
hold on
plot(q1,MSE3);
hold on
plot(q1,MSE4);
hold on
plot(q1,MSE5);
hold on
plot(q1,MSE6);


ylabel('MSE'); 
xlabel('Time'); 
title('MSE with Different Learned R');
legend('Initial R=1','Initial R=0.75','Initial R=0.5','Initial R=0.25','Final Learned R','True R=0.1');
hold off


