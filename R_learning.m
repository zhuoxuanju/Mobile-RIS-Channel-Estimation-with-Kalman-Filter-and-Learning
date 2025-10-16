epoch=60;

I=20;          %ratio between Th and Tg
x1=zeros(1,I);
x2=zeros(1,I);
N=16;           %number of elements
M=8;            %number of antanas
Pt=10;        %transmit power
period=10;
fade_var = 0.1; % variance of  Rayleigh channel

delta=0.01;
eta=0.07;
MSE=zeros(I,period);
MSE2=zeros(I,period);
MSE3=zeros(I,period);
NMSE=zeros(I,period);       %NMSE

sigma_W=sqrt(0.1);
R_real=sigma_W^2*ones(1,epoch+1);
R_learn=zeros(1,epoch+1);
R_learn(1)=1;
R_fixed=R_learn(1)*ones(1,epoch+1);
R_test=zeros(1,epoch);
R_test(1)=R_learn(1)+delta;
sigma_V=0.05;
q=0.001;
h=zeros(N,N,I);
Q=sigma_V^2*eye(N);

phi = sqrtm(Pt) * dftmtx(16)/4;      %get the optimal phi based on the lemma 1
Y=zeros(M,N,I,period);

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
    
    
    
    end
end


for x=1:period
    for b=2:I
    for j=1:N
    W=(1/sqrt(2))*normrnd(0,sigma_W,M,1)+ (1/sqrt(2))*1i * normrnd(0,sigma_W,M,1); % Observation noise
    C(:,:,j)=G*diag(phi(j,:));
    Y1=C(:,:,j)*h(:,j,b)+W;
    Y(:,j,b,x)=Y1;
    end
    end
end

    %h(i,1)=1+1j;
for e=1:epoch


    
for x=1:period           % 5000 Monti Karlo




h_est=zeros(N,1);      %channel estimation

P_pri=zeros(N,N);      %prior covariance
P_post=zeros(N,N);      %posterior covariance
K=zeros(N,M);           %Kalman filter







 

%kalman filter algorithm
P_post=0.1*eye(16);
average_value = mean(h, 'all')/2;
h_est = ones(16,1) * average_value;

MSE(1,x)=(norm(h(:,1,1)-h_est))^2;

for i=2:I
   
    
    for o=1:16
        
       
        Q=q*eye(N)*diag((h_est.*conj(h_est)));
        
        P_pri=P_post+Q;
        K=P_pri*(C(:,:,o))'*inv((C(:,:,o))*P_pri*(C(:,:,o))'+R_learn(e)*eye(M));
        %K=P_pri*(sqrt(Pt)*C_final)'*(1/(sigma_W)^2*eye(128)-1/(sigma_W)^2*(sqrt(Pt)*C_final)*inv((sqrt(Pt)*C_final)'*(sqrt(Pt)*C_final)+sigma_W^2*inv(P_pri))*(sqrt(Pt)*C_final)');
        h_est=h_est+K*(Y(:,o,i,x)-(C(:,:,o))*h_est);
       
        
        P_post=(eye(16)-K*(C(:,:,o)))*P_pri;
    
    end
    
    
    MSE(i,x)=(norm(h(:,16,i)-h_est,2))^2;
    if  e==20
        NMSE=mean(MSE,2);
    end
   



end
end
MSE_mean=mean(MSE,2);
MSE_derive=MSE_mean(I);


for x=1:period           % 5000 Monti Karlo




h_est=zeros(N,1);      %channel estimation

P_pri=zeros(N,N);      %prior covariance
P_post=zeros(N,N);      %posterior covariance
K=zeros(N,M);           %Kalman filter







 

%kalman filter algorithm
P_post=0.1*eye(16);
average_value = mean(h, 'all')/2;
h_est = ones(16,1) * average_value;



MSE2(1,x)=(norm(h(:,1,1)-h_est))^2;

for i=2:I
   
    
    for o=1:16
        
       
        Q=q*eye(N)*diag((h_est.*conj(h_est)));
        
        P_pri=P_post+Q;
        K=P_pri*(C(:,:,o))'*inv((C(:,:,o))*P_pri*(C(:,:,o))'+R_test(e)*eye(M));
        %K=P_pri*(sqrt(Pt)*C_final)'*(1/(sigma_W)^2*eye(128)-1/(sigma_W)^2*(sqrt(Pt)*C_final)*inv((sqrt(Pt)*C_final)'*(sqrt(Pt)*C_final)+sigma_W^2*inv(P_pri))*(sqrt(Pt)*C_final)');
        h_est=h_est+K*(Y(:,o,i,x)-(C(:,:,o))*h_est);
       
        
        P_post=(eye(16)-K*(C(:,:,o)))*P_pri;
    
    end
    
   
    MSE2(i,x)=(norm(h(:,16,i)-h_est,2))^2;
   
   



end
end
MSE2_mean=mean(MSE2,2);
MSE2_derive=MSE2_mean(I);
x1(e)=MSE_derive;


change(e)=(MSE2_derive-MSE_derive)/delta;
R_learn(e+1)=R_learn(e)-eta*change(e);
if R_learn(e+1)<0
    R_learn(e+1)=R_real(e+1)+0.1;
end
R_test(e+1)=R_learn(e+1)+delta;
% if abs(R_learn(e+1)-R_real(e+1))<0.1
%     eta=0.25;
% end

for x=1:period           % 5000 Monti Karlo





h_est=zeros(N,1);      %channel estimation

P_pri=zeros(N,N);      %prior covariance
P_post=zeros(N,N);      %posterior covariance
K=zeros(N,M);           %Kalman filter




% generate Rayleigh fading channel



 

%kalman filter algorithm
P_post=0.1*eye(16);
average_value = mean(h, 'all')/2;
h_est = ones(16,1) * average_value;




MSE3(1,x)=(norm(h(:,1,1)-h_est))^2;

for i=2:I
   
    
    for o=1:16
        
       
        Q=q*eye(N)*diag((h_est.*conj(h_est)));
        
        P_pri=P_post+Q;
        K=P_pri*(C(:,:,o))'*inv((C(:,:,o))*P_pri*(C(:,:,o))'+R_fixed(e)*eye(M));
        %K=P_pri*(sqrt(Pt)*C_final)'*(1/(sigma_W)^2*eye(128)-1/(sigma_W)^2*(sqrt(Pt)*C_final)*inv((sqrt(Pt)*C_final)'*(sqrt(Pt)*C_final)+sigma_W^2*inv(P_pri))*(sqrt(Pt)*C_final)');
        h_est=h_est+K*(Y(:,o,i,x)-(C(:,:,o))*h_est);
       
        
        P_post=(eye(16)-K*(C(:,:,o)))*P_pri;
    
    end
    
    
    MSE3(i,x)=(norm(h(:,16,i)-h_est,2))^2;
    
   



end
end
MSE3_mean=mean(MSE3,2);
MSE3_derive=MSE3_mean(I);
x2(e)=MSE3_derive;



end

period=1000;
for a=1:6
if a==1
    R=1;
elseif a==2
    R=0.75;
elseif a==3
    R=0.5;
elseif a==5
    R=R_learn(epoch);
elseif a==4
    R=0.25;
elseif a==6
    R=0.1;
end



for x=1:period           % 5000 Monti Karlo





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

        h_est=h_est+K*(Y(:,o,i,1)-(C(:,:,o))*h_est);

        P_post=(eye(16)-K*(C(:,:,o)))*P_pri;

    end

    MSE(i,x)=(norm(h(:,16,i)-h_est,2))^2;




end
end
if a==1
    MSE1=mean(MSE,2);
elseif a==2
    MSE2=mean(MSE,2);
elseif a==3
    MSE3=mean(MSE,2);
elseif a==4
    MSE4=mean(MSE,2);
elseif a==5
    MSE5=mean(MSE,2);
elseif a==6
    MSE6=mean(MSE,2);
end
end

q2=1:1:epoch;
q3=1:1:epoch+1;
q1=1:1:I;
plot(q1,NMSE);
pause
plot(q3,R_learn,'r');
ylabel('R'); 
xlabel('Epoch'); 
hold on
plot(q3,R_real,'b');
xlim([0 60]);

title('R with Epoch');
legend('Learning R KF Algorithm','Fixed R KF Algorithm');

hold off
pause

plot(q2,x1,'r');
ylabel('MSE'); 
xlabel('Epoch'); 
hold on
plot(q2,x2,'b');
xlim([0 60]);
title('MSE with Epoch');
hold off

legend('Learning R KF Algorithm','Fixed R KF Algorithm');
pause


plot(q3,R_learn,'r');
ylabel('R'); 
xlabel('Epoch'); 
hold on

xlim([0 60]);
ylim([0 1.2]);
plot(q3,R_fixed,'k');

plot(q3,R_real,'b');
title('R with Epoch');
legend('Learning R KF Algorithm','Fixed R KF Algorithm','True R');

hold off
pause
plot(q1,MSE1);
ylabel('MSE'); 
xlabel('Time'); 
hold on
plot(q1,MSE2);

plot(q1,MSE3);

plot(q1,MSE4);

plot(q1,MSE5);
plot(q1,MSE6);

title('MSE with Different Learned R');
legend('Initial R=1','R=0.75','R=0.5','R=0.25','Final Learned R','True R=0.1');
hold off