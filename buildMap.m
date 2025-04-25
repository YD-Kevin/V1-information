function theta_map=buildMap(rho,N)

    % 生成方向偏好图的二维矩阵
    % 输入:
    %   rho - 密度参数
    %   N   - 网格尺寸
    
    % 输出:
    %   theta_map - N×N的方向偏好矩阵，元素为0到pi的方向
    
    
    
    % 用离散采样近似计算积分的采样参数
    K=360;  %采样数                
    d_alpha=2*pi/K;   %步长     
    alpha=linspace(0, 2*pi, K+1); 
    alpha=alpha(1:end-1); 
    Gamma=2*pi*rand(1,K);  % 随机相位[0, 2π)
    
    % 生成坐标网格（0-based）
    [x1,x2]=meshgrid(0:N-1,0:N-1);
    
    % 计算复数值场phi
    phi=zeros(N, N, 'like', 1i); % 初始化复数矩阵
    for k=1:K
        c=cos(alpha(k));
        s=sin(alpha(k));
        phase=rho*(x1*c+x2*s)+Gamma(k); % 相位计算
        phi=phi+exp(1i*phase)*d_alpha;    % 积分累加
    end
    
    % 计算方向映射
    angle_phi=angle(phi);           % 获取相位角[-pi, pi]
    theta_map=(angle_phi+pi)/2;       %映射到[0,pi]

end