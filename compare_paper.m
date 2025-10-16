I=20;          %ratio between Th and Tg
Th=20e-6;      %coherence time of h

N=16;           %number of elements
M=8;            %number of antanas
Pt=10;        %transmit power
period=1000;
fade_var = 0.1; % variance of  Rayleigh channel


MSE=zeros(I,period);
MSE2=zeros(I,period);
MSE3=zeros(I,period);


%sigma_W=0.5;

sigma_V=0.05;
q=0.001;
h=zeros(N,N,I);
Q=sigma_V^2*eye(N);

phi = sqrtm(Pt) * dftmtx(16)/4;      %get the optimal phi based on the lemma 1
Y=zeros(M,N,I);
Y2=zeros(N,M,I);
D=zeros(N,M,I);
v=300/3.6;        %speed
f=2600000000;     %carrier frequency
f_dopshift=v*f/300000000;     %maximum doppler shift frequency
alpha=besselj(0,2*pi*f_dopshift*Th);      %tried to use zero-order Bessel function of the first kind
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

for e=1:2
if e==1
    sigma_W=sqrt(0.5);
elseif e==2
    sigma_W=sqrt(0.1);
end
R=sigma_W^2;
for b=2:I
    for j=1:N
    n=(1/sqrt(2))*normrnd(0,sigma_V,N,1)+ (1/sqrt(2))*1i * normrnd(0,sigma_V,N,1); % Inside noise
    
    h1=h1+h1.*n;
    
    h(:,j,b)=h1;
    
     W=(1/sqrt(2))*normrnd(0,sigma_W,M,1)+ (1/sqrt(2))*1i * normrnd(0,sigma_W,M,1); % Observation noise
     C(:,:,j)=G*diag(phi(j,:));
     Y1=C(:,:,j)*h(:,j,b)+W;
     Y(:,j,b)=Y1;
     W2(:,j,b)=W;
    end
    D(:,:,b)=diag(h(:,N,b))*G';
    Y2(:,:,b)=phi.'*D(:,:,b)+W2(:,:,b).';
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
h_est = ones(16,1) * average_value;


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


end



for x=1:period
h_est2=zeros(N,1);      %channel estimation

P_pri=zeros(N,N);      %prior covariance
P_post=zeros(N,N);      %posterior covariance
K=zeros(N,N);           %Kalman filter




% generate Rayleigh fading channel



 

%kalman filter algorithm
P_post=M*eye(16);
average_value = mean(h(:,1,1), 'all')/2;
h_est2 = zeros(16,1) * average_value;
D_est=zeros(N,M);

MMSE(1,x)=(norm(h(:,1,1)-h_est2))^2;

for i=2:I

D_est=alpha*D_est;
P_pri=alpha^2*P_post+M*(1-alpha^2)*eye(N);
K=P_pri*phi'*inv((phi.'*P_pri*phi'+M*eye(N)));
D_est=D_est+K*(Y2(:,:,b)-phi.'*D_est);
P_post=(eye(N)-K*phi')*P_pri;
MMSE(i,x)=(norm(D(:,:,b)-D_est,'fro'))^2;
end
end

if e==1
    MMSE1=mean(MMSE,2);
elseif e==2
    MMSE2=mean(MMSE,2);


end


end


q1=1:1:I;
plot(q1,MSE1);
hold on
plot(q1,MSE2);
hold on
plot(q1,MMSE1);
hold on
plot(q1,MMSE2);



ylabel('MSE'); 
xlabel('Time'); 
title('MSE with Time');
legend('Algorithm 1, R=0.5','Algorithm 1, R=0.1','[15], R=0.5','[15], R=0.1');
hold off



