function gabor = generateGabor(lambda, theta)
    %生成Gabor滤波器
    % 输入：
    % lambda-波长
    % theta-方向角度

    % 输出： 滤波核矩阵gabor
    
    lambda = round(lambda); % 确保边长为整数
    if lambda <= 0
        error('lambda 必须为正整数');
    end
    sigma = sqrt(lambda); % 高斯项标准差（方差为 lambda）
    theta_rad = deg2rad(theta); % 角度转弧度
    
    % --- 生成坐标网格 ---
    s = 2 * lambda + 1; % 滤波器边长
    [X, Y] = meshgrid(-lambda:lambda, -lambda:lambda);
    
    % --- 旋转坐标（方向变换）---
    x_theta =  X * cos(theta_rad) + Y * sin(theta_rad);
    y_theta = -X * sin(theta_rad) + Y * cos(theta_rad);
    
    % --- 计算 Gabor 函数核心 ---
    gamma = 0.5; % 椭圆率（默认值）
    psi = 0; % 相位偏移（默认值）
    
    % 高斯调制项
    gaussian = exp(-(x_theta.^2 + (gamma^2) * y_theta.^2) / (2 * sigma^2));
    
    % 余弦载波项
    cosine = cos(2 * pi * x_theta / lambda + psi);
    
    % Gabor 滤波器
    gabor = gaussian .* cosine;
    
    % --- 归一化到 [-1, 1] ---
    gabor = gabor / max(abs(gabor(:))); % 确保最大绝对值为 1
end