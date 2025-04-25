function makeGif(list,time)

time=min(time,size(squeeze(list(:,1,1)))); %帧数量
filename = 'animation.gif';  % 输出文件名
delay_time = 0.1;            % 每帧延迟时间（秒）
loop_count = Inf;            


for k = 1:time
    % 提取当前帧并归一化到 [0, 1]（假设原始数据在合理范围内）
    current_frame = squeeze(list(k, :, :));
    normalized_frame = mat2gray(current_frame); 
    
    
    uint8_frame = im2uint8(normalized_frame);% 转换为 uint8 格式（0-255）
    
    
    cmap = colormap(gray(256));% 生成颜色映射（灰度）


    % 写入GIF文件
    if k == 1
        imwrite(uint8_frame, cmap, filename, 'gif', ...
            'DelayTime', delay_time, 'Loopcount', loop_count);
    else
        imwrite(uint8_frame, cmap, filename, 'gif', ...
            'DelayTime', delay_time, 'WriteMode', 'append');
    end
end

disp(['生成GIF，已保存为: ', filename]);
end