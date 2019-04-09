clear; close all;
[y,Fs] = audioread('NMGm.wav');
[yn,Fsn] = audioread('noise.wav');
[yc,Fsc] = audioread('o1.wav');




fig1 = figure(1);
set(fig1,'Position',[0,0,550,250]);
stem(y,'Marker','none');
grid on;
xlabel('$$k$$','Interpreter', 'Latex','FontSize',16)
ylabel('$$x(k)$$','Interpreter', 'Latex','FontSize',16)

fig2 = figure(2);
set(fig2,'Position',[0,0,550,250]);
stem(yn,'Marker','none');
grid on;
xlabel('$$k$$','Interpreter', 'Latex','FontSize',16)
ylabel('$$x(k)$$','Interpreter', 'Latex','FontSize',16)

fig3 = figure(3);
set(fig3,'Position',[0,0,550,250]);
stem(yc,'Marker','none');
grid on;
xlabel('$$k$$','Interpreter', 'Latex','FontSize',16)
ylabel('$$x(k)$$','Interpreter', 'Latex','FontSize',16)

num = 8000;

fig4 = figure(4);
c = fft(y,num);
t = (0:1/num:0.5-1/num)';
set(fig4,'Position',[0,0,550,250]);
stem(t,abs(c(1:num/2)),'Marker','none');
grid on;
xlabel('$$f$$','Interpreter', 'Latex','FontSize',16)
ylabel('$$X(f)$$','Interpreter', 'Latex','FontSize',16)

fig5 = figure(5);
cn = fft(yn,num);
tn = (0:1/num:0.5-1/num)';
set(fig5,'Position',[0,0,550,250]);
stem(tn,abs(cn(1:num/2)),'Marker','none');
grid on;
xlabel('$$f$$','Interpreter', 'Latex','FontSize',16)
ylabel('$$X(f)$$','Interpreter', 'Latex','FontSize',16)

fig6 = figure(6);
cc = fft(yc,num);
tc = (0:1/num:0.5-1/num)';
set(fig6,'Position',[0,0,550,250]);
stem(tc,abs(cc(1:num/2)),'Marker','none');
grid on;
xlabel('$$f$$','Interpreter', 'Latex','FontSize',16)
ylabel('$$X(f)$$','Interpreter', 'Latex','FontSize',16)

fig7 = figure(7);
set(fig7,'Position',[0,0,550,250]);
stem(t,abs(cn(1:num/2))-abs(c(1:num/2)),'Marker','none');
grid on;
xlabel('$$f$$','Interpreter', 'Latex','FontSize',16)
ylabel('$$X(f)$$','Interpreter', 'Latex','FontSize',16)

fileid = fopen('filter_data.dat');
num = fread(fileid, 1, 'uint32');
a(:,1) = fread(fileid, num, 'double');

saveas(fig1,'figdata/audio_x','pdf');
saveas(fig2,'figdata/audio_xn','pdf');
saveas(fig3,'figdata/audio_xc','pdf');
saveas(fig4,'figdata/audio_dft_x','pdf');
saveas(fig5,'figdata/audio_dft_xn','pdf');
saveas(fig6,'figdata/audio_dft_xc','pdf');
saveas(fig7,'figdata/audio_dft_xd','pdf');
% 
% 
% num = 8000;
% % plot(y)e
% close all;
% figure(1)
% c = fft(y,num);
% t = (0:1/num:1-1/num)';
% stem(t,abs(c),'Marker','none');
% 

% figure(2)
% [y_r,Fs_r] = audioread('NMGm.wav');
% c_r = fft(y_r,num);
% stem(t,abs(c_r),'Marker','none');
% 
% figure(3)
% stem(t,abs(c)-abs(c_r),'Marker','none');
% grid on;
% 
% 
% 
% figure(4)
% [y_nc,Fs_nc] = audioread('o1.wav');
% c_nc = fft(y_nc,num);
% stem(t,abs(c_nc),'Marker','none');
% 
% figure(5)
% subplot(3,1,1)
% plot(y_r)
% subplot(3,1,2)
% plot(y)
% subplot(3,1,3)
% plot(y_nc)
% 
% 
% 
