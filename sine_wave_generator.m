%% 生成正弦信号
function y=sine_wave_generator(f)
global DataLong;
A=1;%幅度
w=2*pi*f;           
fs=8000;%采样频率
dt=1/fs;%采样间隔
t=0:dt:dt*(DataLong-1);
y=A*sin(w*t);%正弦信号
end