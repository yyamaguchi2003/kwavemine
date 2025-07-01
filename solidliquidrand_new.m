clearvars;
close all;
DATA_CAST='gpuArray-single';
config = jsondecode(fileread('/home/yamaguchi/document/tutorial/kwavemine/config2d.json'));
Nx = 4*config.grid.Nx;
Ny = 4*config.grid.Ny;
dx = 0.25*config.grid.dx;
dy = 0.25*config.grid.dy;
kgrid = kWaveGrid(Nx, dx, Ny, dy);
save_path = config.save_path;
t_end = config.simulation.t_end
CFL = config.simulation.CFL

% -------------------------------------------------------------------------
% 3) 媒質パラメータ
% -------------------------------------------------------------------------
% 水のパラメータ
medium.sound_speed = config.medium.water.sound_speed;
medium.density     = config.medium.water.density;
medium.alpha_coeff = config.medium.water.alpha_coeff;
medium.alpha_power = config.medium.water.alpha_power;
medium.BonA = config.medium.water.BonA;

% ガラスのパラメータ
vinyl.sound_speed = config.medium.vinyl.sound_speed;
vinyl.density     = config.medium.vinyl.density;
vinyl.alpha_coeff = config.medium.vinyl.alpha_coeff;

%空気のパラメータ
air.sound_speed = config.medium.air.sound_speed;
air.density     = config.medium.air.density;
air.alpha_coeff = config.medium.air.alpha_coeff;
air.BonA        = config.medium.air.BonA;

%ガラスビーズのパラメータ
glass.sound_speed = config.medium.glass.sound_speed;
glass.density = config.medium.glass.density;
glass.alpha_coeff = config.medium.glass.alpha_coeff;
glass.BonA = config.medium.glass.BonA;



distance_pipe_source = config.source.distance_pipe_source;
source_diameter = config.source.diameter;
outer_radius = config.pipe.outer_radius;  % 外側の半径
inner_radius = config.pipe.inner_radius; % 内側の半径

num_particle = config.fraction.number_particle;
fraction_particle = config.fraction.fraction_particle;
diameter_particle = config.fraction.diameter_particle;
pipe_height = num_particle*4*diameter_particle^3/(6*inner_radius^2*fraction_particle);
pipe_height = pipe_height - diameter_particle;
gap_particles = config.fraction.gap_particles;

gap_pipe = config.fraction.gap_pipe_particle;   %パイプと粒子の間の距離
d_max = inner_radius - diameter_particle/2;

% -------------------------------------------------------------------------
% 5) 塩ビパイプの設定
% -------------------------------------------------------------------------
% パイプの中心位置
center_x = Nx/2 + config.pipe.center_x_offset;
center_y = Ny/2 + config.pipe.center_y_offset;

% 媒質パラメータのマスクを作成
medium.sound_speed = medium.sound_speed * ones(Nx, Ny);
medium.density = medium.density * ones(Nx, Ny);
medium.alpha_coeff = medium.alpha_coeff * ones(Nx, Ny);
medium.BonA = medium.BonA * ones(Nx, Ny);

% 塩ビパイプのパラメータを設定
for i = 1:Nx
    for j = 1:Ny
        if dx*sqrt((i-Nx/2)^2+(j-Ny/2)^2)<= outer_radius & dx*sqrt((i-Nx/2)^2+(j-Ny/2)^2)>= inner_radius
            medium.sound_speed(i, j) = vinyl.sound_speed;
            medium.density(i, j) = vinyl.density;
            medium.alpha_coeff(i, j) = vinyl.alpha_coeff;
            medium.BonA(i, j) = 0;
        end
    end
end

% ガラスビーズのパラメータ

%fpsolid = fopen(fullfile(save_path, 'solidliquidrandposition.txt'), 'r');
%center_a = textscan(fpsolid, '%f %f %f');
%fclose(fpsolid);

%center_all = [center_a{1}*1e-3, center_a{2}*1e-3, center_a{3}*1e-3];

% 新しく乱数を生成する場合↓

center_all = zeros(num_particle,3);
rng('shuffle');

for i=1:num_particle
    while 1
        x_rand = sqrt(-2*log(rand()));
        if x_rand >= 2
            continue;
        end
        y_rand = 2*pi*rand();
        z_pole = pipe_height*rand()-0.5*pipe_height;
        center_all(i,:) = [d_max/2*x_rand*cos(y_rand) d_max/2*x_rand*sin(y_rand) z_pole];
        gap_min = inf;
        distance_pars= 0;
        if i>1
            for j=1:i-1
                distance_pars = sqrt((center_all(i,1)-center_all(j,1))^2+(center_all(i,2)-center_all(j,2))^2+(center_all(i,3)-center_all(j,3))^2);
                if gap_min > distance_pars
                    gap_min = distance_pars;
                end
            end
        end
        if gap_min> diameter_particle
            disp(center_all(i,:));
            break;
        end
    end
end

% x y z r^2 の形式でtxtに保存

fp = fopen(fullfile(save_path, 'solidliquidrand4gnuplot2.txt'), 'w');
for i=1:num_particle
    if diameter_particle^2/4 > center_all(i,3)^2 
        fprintf(fp, "set object %d circle at %.3f, %.3f size %.3f lc 'black' lw 2\n", i+2, center_all(i,2)*1e3, 11.8+center_all(i,1)*1e3, sqrt(diameter_particle^2/4-center_all(i,3)^2)*1e3);
    end
end
fprintf(fp, "No 11.8 offset from here\n")
for i=1:num_particle
    if diameter_particle^2/4 > center_all(i,3)^2 
        fprintf(fp, "set object %d circle at %.3f, %.3f size %.3f lc 'black' lw 2\n", i+2, center_all(i,2)*1e3, center_all(i,1)*1e3, sqrt(diameter_particle^2/4-center_all(i,3)^2)*1e3);
    end
end
fclose(fp);

fp2 = fopen(fullfile(save_path, 'solidliquidrandposition2.txt'), 'w');
for i=1:num_particle
    if diameter_particle^2/4 > center_all(i,3)^2 
        fprintf(fp, "%.3f %.3f %.3f\n", center_all(i,1)*1e3, center_all(i,2)*1e3, center_all(i,3)*1e3);
    end
end
fclose(fp2);

for i=1:Nx
    for j=1:Ny
        for k=1:num_particle
            if sqrt(((i-Nx/2)*dx-center_all(k,1))^2+((j-Ny/2)*dx-center_all(k,2))^2+center_all(k,3)^2) < diameter_particle/2
                medium.sound_speed(i, j) = glass.sound_speed;
                medium.density(i, j) = glass.density;
                medium.alpha_coeff(i, j) = glass.alpha_coeff;
                medium.BonA(i, j) = 0;
            end
        end
    end
end

% -------------------------------------------------------------------------
% 6) シミュレーション時間配列の作成
% -------------------------------------------------------------------------
kgrid.makeTime(medium.sound_speed, CFL, t_end);
distance_pipe_source = distance_pipe_source/dy;
source_diameter = source_diameter/dx;

source.p_mask = zeros(Nx, Ny);
source.p_mask(round(Nx/2-source_diameter/2):round(Nx/2+source_diameter/2), round(Ny/2-distance_pipe_source))= 1;
source_signal = zeros(size(kgrid.t_array));
source_frequency = config.source.tone_burst_freq;
T_prf = 1 / config.source.prf;
t_array = kgrid.t_array;

for n = 0:config.source.max_n
    t_start = n * T_prf;
    t_end = t_start + config.source.pulse_length;
    
    if t_start > kgrid.t_array(end)
        break;
    end
    
    idx_on = (t_array >= t_start) & (t_array < t_end);
    source_signal(idx_on) = config.source.magnitude * (sin(2*pi * source_frequency * (t_array(idx_on) - t_start)));
%source_signal = toneBurst(1.0/kgrid.dt, config.source.frequency, 4, 'Envelope', 'Rectangular');
end
source.p = source_signal;

% -------------------------------------------------------------------------
% 8) センサーの設定
% -------------------------------------------------------------------------
sensor.mask = zeros(Nx, Ny);
sensor_x = Nx/2;
sensor_diameter = config.sensor.diameter;
sensor_diameter = sensor_diameter/dx;
sensor_y_offset = config.sensor.y_offset/dy;
sensor_y = round(Ny/2 + sensor_y_offset);
sensor.mask(round(Nx/2-sensor_diameter/2):round(Nx/2+sensor_diameter/2), sensor_y) = 1;
sensor_y = round(Ny/2 - sensor_y_offset);
sensor.mask(round(Nx/2-sensor_diameter/2):round(Nx/2+sensor_diameter/2), sensor_y) = 1;
sensor.record = {'p'};

% -------------------------------------------------------------------------
% 9) シミュレーションのオプション設定
% -------------------------------------------------------------------------
input_args = {
    'PMLInside', true, 'PlotPML', true, ...
    'PMLSize', config.simulation.pml_size, ...
    'RecordMovie', false, ...
    'PlotFreq', 50, ...
    'MovieName', fullfile(save_path, 'vinylpipe_exp.avi'), ...
    'DataCast', DATA_CAST, ...
    };

% -------------------------------------------------------------------------
% 10) シミュレーション実行
% -------------------------------------------------------------------------
sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});

% -------------------------------------------------------------------------
% 11) 結果の可視化
% -------------------------------------------------------------------------
len = length(sensor_data.p(:,1));
ref = sensor_data.p(1:len/2,:);
trans = sensor_data.p(len/2+1:len,:);
figure;
plot(kgrid.t_array*1e6, mean(ref)*1e-3);
set(gca, 'fontsize', 15); %文字サイズ変更
set(gca, 'XTick', [0:0.05e3:0.2e3]);
set(gca, 'YTick', [-150:150:150]);
set(gca, 'ylim', [-150, 150]);
xlabel('Time [μs]');
ylabel('Pressure [kPa]');
title('Pressure at TDX1 with a Glass Bead');
saveas(gcf, fullfile(save_path, 'solidliquidrand_ref2.png')); 

ref = mean(ref);

fileid = fopen(fullfile(save_path,"solidliquidrand2_TDX1.txt"), 'w');
for t = 1:round(length(kgrid.t_array));
    fprintf(fileid, "%.5e %.5e\n", kgrid.t_array(t), ref(1,t));
end
fclose(fileid);

figure;
plot(kgrid.t_array*1e6, mean(trans)*1e-3);
set(gca, 'fontsize', 15); %文字サイズ変更
set(gca, 'XTick', [0:0.05e3:0.2e3]);
set(gca, 'YTick', [-150:150:150]);
set(gca, 'ylim', [-150, 150]);
xlabel('Time [μs]');
ylabel('Pressure [kPa]');
title('Pressure at TDX3 with a Glass Bead');
saveas(gcf, fullfile(save_path, 'solidliquidrand_trans2.png')); 

trans = mean(trans);

fileid2 = fopen(fullfile(save_path,"solidliquidrand2_TDX3.txt"), 'w');
for t = 1:round(length(kgrid.t_array));
    fprintf(fileid2, "%.5e %.5e\n", kgrid.t_array(t), trans(1,t));
end
fclose(fileid2);

% kspaceFirstOrder2D実行後にGPU配列をCPUに集約
sensor_data_cpu = structfun(@gather, sensor_data, 'UniformOutput', false);

% v7.3 形式で.matファイル保存
save(fullfile(save_path, 'solidliquidrand2.mat'), 'sensor_data_cpu', '-v7.3');