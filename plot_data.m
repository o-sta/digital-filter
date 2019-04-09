clear; close all;

fileid = fopen('filter_BEF.dat');
num = fread(fileid, 1, 'uint32');
a(:,1) = fread(fileid, num, 'double');
%窓関数法　フィルタ係数のプロット
fig1=figure(1);
stem(a(:,1),'Marker','none','LineWidth',1)
xlabel('$$n$$','Interpreter', 'Latex','FontSize',16)
ylabel('$$h_n$$','Interpreter', 'Latex','FontSize',16)
grid on;

%窓関数法　dBゲインとゲインのプロット
point_num = 10000;
l = pi/point_num;
data_matlab = zeros(point_num+1,2);
for point = 0:1:point_num
    res = 0.0;
    n = 0.0;
    for v = (num-1)/2+1:1:num
        if (n==0)
            res = res + a(v,1)*cos(n*(point*l));
        else
            res = res + 2*a(v,1)*cos(n*(point*l));
        end
        n = n + 1;
    end
    data_matlab(point+1,:) = [point*l res];
end

% 
% fig2=figure(2);
% fileid = fopen('gain_data.dat');
% num = fread(fileid, 1, 'uint32');
% data = fread(fileid, [3,num], 'double');
% plot(data(1,:),data(2,:),'LineWidth',1);
% hold on;
% % plot(data_matlab(:,1),data_matlab(:,2),'LineWidth',1);
% grid on;
% xlabel('$$\omega$$','Interpreter', 'Latex','FontSize',16)
% ylabel('$$H(\omega)$$','Interpreter', 'Latex','FontSize',16)

fig3=figure(3);
% plot(data(1,:),data(3,:),'LineWidth',1);
hold on;
plot(data_matlab(:,1),20*log10(abs(data_matlab(:,2))),'LineWidth',1);
grid on;
xlabel('$$\omega$$','Interpreter', 'Latex','FontSize',16)
ylabel('$$20\log{|H(\omega)|}$$','Interpreter', 'Latex','FontSize',16)

% 
% saveas(fig1,'figdata/h2','pdf');
% saveas(fig2,'figdata/gain2','pdf');
% saveas(fig3,'figdata/db2','pdf');




