clc;
close all;
clear all;
%n=input('enter the satnumber');
nn=20;
load gpsinp_730ms;
a(1:3*5714)=gpsinp_730ms(1:3*5714);
% a(1:20*5714)=int_data(2484:20*5714+2483);
%fs=input('enter the sample frequency');
fs=5714000;
 %Ts=1/5800000;
%n_samples=input('enter number of samples');
n_samples=3*5714;
%chiprate=input('enter the chiprate');
chiprate=1023000;
%rif=input('enter the rf if');
rif=4309000;
%t_sec=input('enter the number of seconds');
% t_sec=.01;
%fc=input('enter the carrieer frequency');
fc=1405000;
% num_sec=input('enter the number of seconds');
 %sat_id=input('enter the satellite number');
N=200;M = 30 ;
for k = 1:86;
X(1:N)=gpsinp_730ms((k-1)*N+1:(k)*N);
t = (1:N)';

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

RC=zeros(N,M);
for m=1:M
  buf=PC(:,m)*RHO(:,m)'; 
  buf=buf(end:-1:1,:);
  for n=1:N 
    RC(n,m)=mean( diag(buf,-(N-M+1)+n) );  
  end
end;
D((k-1)*N+1:(k)*N)=sum(RC(:,1:2),2);
end
E(1:3*5714)=D(1:3*5714);
phase=[2 6;3 7;4 8;5 9;1 9;2 10;1 8;2 9;3 10;2 3;3 4;5 6;6 7;7 8;8 9;9 10;1 4;2 5;3 6;4 7;5 8;6 9;1 3;4 6;5 7;6 8;7 9;8 10;1 6;2 7;3 8;4 9];
Ts=1/fs;
g1=[-1 -1 -1 -1 -1 -1 -1 -1 -1 -1];
g2=g1;
s1=phase(nn,1);
s2=phase(nn,2);
tmp=0;
for i=1:1023;
    g(i)=g2(s1)*g2(s2)*g1(10);
    tmp=g1(1);
    g1(1)=g1(3)*g1(10);
    g1(2:10)=[tmp g1(2:9)];
    tmp=g2(1);
    g2(1)=g2(2)*g2(3)*g2(6)*g2(8)*g2(9)*g2(10);
    g2(2:10)=[tmp g2(2:9)];
  end;  
  l=1:n_samples; 
 ca1(l)=g(mod(ceil(l*chiprate/fs)-1,1023)+1);
 
 ca2 = - ca1;
 fft_ca = fft(ca1);
 fft_ca2 = fft(ca2);
%  for i = 1: 1: n_samples;
%      
%      fft2_ca(i) = fft_ca(n_samples-i+1);
%      
%  end;
%  
% ca2 = ifft(fft2_ca);
%  plot(abs(fft_ca));
%  
%  plot(abs(fft2_ca));
 
%  plot(abs(ifft(fft_ca .* fft_ca)));
 
 
  
 i=1071;

   ca2=[ca1((1+(i-1)*5):3*5714) ca1(1:(i-1)*5)];
% for ab = 2:1:15;
%       ca2(1:n_samples)=[ca2(2:n_samples) ca2(1)];
phase_count=0;
cenfreq_num=1405470*2*pi/5714000;
loop_output_num=0;
  for i=1:1:3*5714;
    phase_count=phase_count+cenfreq_num+loop_output_num;
    sin1(i)=sin(phase_count);
    cos1(i)=cos(phase_count);
     if (sin1(i)>0)
         if (sin1(i)>0.5)
            sin2(i)=2;
         else
             sin2(i)=1;
         end;
     else
         if (sin1(i)<-0.5)
            sin2(i)=-2;
         else
             sin2(i)=-1;
         end;
     end;
     if (cos1(i)>0)
         if (cos1(i)>0.5)
            cos2(i)=2;
         else
             cos2(i)=1;
         end;
     else
         if (cos1(i)<-0.5)
            cos2(i)=-2;
         else
             cos2(i)=-1;
         end;
     end;
    
  end;
  sin3=ca2.*sin2;
  cos3=ca2.*cos2;
  for j=1:256;
    data1(j)=sum(E((1+(j-1)*50):j*50).*sin3((1+(j-1)*50):j*50));
    data2(j)=sum(E((1+(j-1)*50):j*50).*cos3((1+(j-1)*50):j*50));
%      disp(j)
end;
 z=complex(data1,data2);
 outfft=fft(z);
 c=abs(outfft);
%  abc(ab)= max(c);
% end;
 plot(c);
 xlabel('doppler frequency');
 ylabel('amplitude');