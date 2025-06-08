clc; 
close all;
clear all;
n_samples=10*5714;
fs=5714000;
load gpsinp_730ms
a(1:10*5714)=gpsinp_730ms(1:10*5714);


%  n=9;
%  n=24;
% n=16;
n=20;
% n=23;
phase=[2 6;3 7;4 8;5 9;1 9;2 10;1 8;2 9;3 10;2 3;3 4;5 6;6 7;7 8;8 9;9 10;1 4;2 5;3 6;4 7;5 8;6 9;1 3;4 6;5 7;6 8;7 9;8 10;1 6;2 7;3 8;4 9];
s1=phase(n,1);
s2=phase(n,2);
tmp=0;
chiprate=1023000;
g1=[-1 -1 -1 -1 -1 -1 -1 -1 -1 -1];
g2=g1;
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
       for i=1:1142;
%   i=836;
% i=973;
% i=1071;
%    i=1071;
%  i=647;
 ca2=[ca1((1+(i-1)*5):10*5714) ca1(1:(i-1)*5)];
%      ca2=[ca2(2:1*5714) ca2(1:1)];
%         ca2=[ca2((1*5714)-1:1*5714) ca2(1:(1*5714)-2)];
realdata=ca2.*a;
z=fft(realdata);
  c=abs(z);
  max_corr(i)=(max(c));
%         f2=5714000*(1:5*5714)/(10*5714);
%   subplot(3,1,1);
%        plot(f2,c(1,(1:5*5714)));
%        xlabel('frequency');
%        ylabel('amplitude');
%  plot(c);
% plot(c);
        end;
disp('maximum=');
 temp1=0;
 for i=1:1:1142;
     if temp1<max_corr(i);
         temp1=max_corr(i);
         i
     end ;
     end;
% sat_max(n)=max(max_corr)
%    end;
% for n=1:32
%     disp('satid=');
%     disp(n);
%     disp(sat_max(n));
% end;
%plot(abs(z));