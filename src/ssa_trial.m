%%function (D) = ssa_trail(X)

clear;

M = 30;   
N = 200;   
%T = 22;    
%stdnoise = 1; 

load gpsinp_730ms
X(1:10*5714)=gpsinp_730ms(1:10*5714);
t = (1:N)';
            

figure(1);
set(gcf,'name','Time series X');
clf;
plot(t,X,'b-');



covX = xcorr(X,M-1,'unbiased');
Ctoep=toeplitz(covX(M:end));

C=Ctoep;

[RHO,LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);              
[LAMBDA,ind]=sort(LAMBDA,'descend');
RHO = RHO(:,ind);                   
Y=zeros(N-M+1,M);
for m=1:M    
  Y(:,m) = X((1:N-M+1)+m-1);
end;
PC = Y*RHO;

figure(4);
set(gcf,'name','Principal components PCs')
clf;
for m=1:4
  subplot(4,1,m);
  plot(t(1:N-M+1),PC(:,m),'k-');
  ylabel(sprintf('PC %d',m));
  ylim([-10 10]);
end;

RC=zeros(N,M);
for m=1:M
  buf=PC(:,m)*RHO(:,m)'; 
  buf=buf(end:-1:1,:);
  for n=1:N 
    RC(n,m)=mean( diag(buf,-(N-M+1)+n) );  
  end
end;

figure(5);
set(gcf,'name','Reconstructed components RCs')
clf;
for m=1:4
  subplot(4,1,m);
  plot(t,RC(:,m),'r-');
  ylabel(sprintf('RC %d',m));
  ylim([-1 1]);
end;

figure(6);
set(gcf,'name','Original time series X and reconstruction RC')
clf;
subplot(2,1,1)
plot(t,X,'b-',t,sum(RC(:,:),2),'r-');
legend('Original','Complete reconstruction');

subplot(2,1,2)
plot(t,X,'b','LineWidth',2);
D = sum(RC(:,1:2),2);
D = D';
plot(t,X,'b-',t,sum(RC(:,1:2),2),'r-');
legend('Original','Reconstruction with RCs 1-2');
end