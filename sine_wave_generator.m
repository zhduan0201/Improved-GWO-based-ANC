%% ���������ź�
function y=sine_wave_generator(f)
global DataLong;
A=1;%����
w=2*pi*f;           
fs=8000;%����Ƶ��
dt=1/fs;%�������
t=0:dt:dt*(DataLong-1);
y=A*sin(w*t);%�����ź�
end