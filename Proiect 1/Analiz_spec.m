function [P1, f] = Analiz_spec(x,Fs,ti)
y=fft(x);
L=length(y);
P2 = y/L;
P1 = 2*P2(1:fix(L/2)+1);
f = Fs*(0:fix(L/2))/L;
figure;
plot (f,abs(P1));
title(ti)
xlabel('f (Hz)')
ylabel('|P1(f)|')
end