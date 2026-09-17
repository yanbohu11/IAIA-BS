function [fitness,bX] = fobj(X,opt)
[~,dim]=size(X);
bX=zeros(1,dim);
                                                                                          
%二进制映射
for j=1:dim
    x=1./(1+exp(-X/2));
    bX(j)=double(x(j)<rand());
end

q=opt.band_num;
N=length(bX);
indices = ceil((1:N) * q / N);
groups = cell(1, q); 
extracted_values = cell(1, q);
selected_idx = zeros(1, 1); 
for i = 1:q
  groups{i} = bX(indices == i);  % 按组索引提取元素 
end
% group_lengths = cellfun(@length, groups);
% cumulative_lengths = [0, cumsum(group_lengths)];
for i=1:q
    count=sum(groups{i});
    if count==1
       
    elseif count==0
        extracted_values{i} = opt.LD(indices==i);
        [~, max_idx] = max(extracted_values{i});
        groups{i}(max_idx)=1;
    else
        indices2 = find(groups{i} == 1);   %元素值为1的索引
        candidates = indices2(randperm(length(indices2), 2));
        z1=opt.LD(candidates(1));
        z2=opt.LD(candidates(2));
        values = [z1, z2];
        % [rows, cols] = meshgrid(indices2, indices2);
        % selected_elements = opt.matrix(sub2ind(size(opt.matrix), rows(:), cols(:)));
        % indices2lenth = length(indices2);    %索引长度
        % reshaped_elements = reshape(selected_elements, indices2lenth, indices2lenth);
        % column_sums = sum(reshaped_elements, 1);
        % candidates = indices2(randperm(length(indices2), 2));  %candidates =  randsample(indices2, 2);                        
        % rank1 = indices2 == candidates(1);
        % rank2 = indices2 == candidates(2);
        % positions = arrayfun(@(t) find(indices2== t),  candidates);          % position =find(indices2 == candidates);
        % values = [column_sums(rank1), column_sums(rank2)];
        [~, best] = max(values);
        selected_idx = candidates(best);
        
        new_x = zeros(size(groups{i}));
        new_x(selected_idx) = groups{i}(selected_idx);
        groups{i}=new_x;
    end 
    
end
bX = zeros(1, N);
for i = 1:q
    bX(indices == i) = groups{i};
    % disp(bX);
end


% [cor,gt] = get_dataset_name(2);
% [H, W, ~] = size(cor);
% labels=gt;
% selected_bands=find(X);
% nB=length(selected_bands);
% cors=cor(:,:,selected_bands);
% cors1 = reshape(cors, [], nB);
% labels = labels(:); 
% classes = unique(labels);
% classes(classes == 0) = []; % 0 可能是无效标签
% num_classes = length(classes);
% 
% mu = arrayfun(@(c) mean(cors1(labels == c, :), 1), classes, 'UniformOutput', false);
% sigma = arrayfun(@(c) var(cors1(labels == c, :), 0, 1), classes, 'UniformOutput', false);
% mu = cat(1, mu{:});  % 变成 num_classes × bands 矩阵
% sigma = cat(1, sigma{:}); % 变成 num_classes × bands 矩阵
% 
% %4. 计算 JM 距离矩阵（无 for 循环）
% i_idx = repmat((1:num_classes)', 1, num_classes); % 行索引
% j_idx = i_idx'; % 列索引
% mu_diff = mu(i_idx, :) - mu(j_idx, :); % 均值差值 (num_classes × num_classes × bands)
% sigma_avg = 0.5 * (sigma(i_idx, :) + sigma(j_idx, :)); % 平均方差
% 
% B = 0.125 * (mu_diff.^2 ./ sigma_avg) + 0.5 * log(sigma_avg ./ sqrt(sigma(i_idx, :) .* sigma(j_idx, :)));
% JM_matrix = 2 * (1 - exp(-sum(B, 3))); % 计算 JM 距离矩阵
% JMtotal = sum(JM_matrix(:)); 




% selected_LD=bX.*opt.LD;             %%f1       
% [~,sort_indices]=sort(selected_LD);
% sort_indices=flip(sort_indices);  %反转元素顺序
% sort_indices=sort_indices(1:opt.band_num);  %sort_indices为前30个的下标
% temp_array=zeros(1,dim);   %[1~103]
% temp_array(sort_indices)=1;  %选前30个设为1
% bX=bX.*temp_array;      %%f1 不改变bX
%disp(num2str(sum(bX)));



% [H, W, ~] = size(opt.cor);
% labels=opt.gt;
% selected_bands=find(X);
% nB=length(selected_bands);
% cors=opt.cor(:,:,selected_bands);
% cors1 = reshape(cors, [], nB);
% 
% mu = zeros(opt.num_classes, nB);
% sigma = zeros(opt.num_classes, nB);
% for c = 1:opt.num_classes
%     class_idx = labels == opt.classes(c);
%     mu(c, :) = mean(cors1(class_idx, :), 1);
%     sigma(c, :) = var(cors1(class_idx, :), 0, 1);
% end
% 
% 
% % mu = arrayfun(@(c) mean(cors1(labels == c, :), 1), classes, 'UniformOutput', false);
% % sigma = arrayfun(@(c) var(cors1(labels == c, :), 0, 1), classes, 'UniformOutput', false);
% % mu = cat(1, mu{:});  % 变成 num_classes × bands 矩阵
% % sigma = cat(1, sigma{:}); % 变成 num_classes × bands 矩阵
% 
% %4. 计算 JM 距离矩阵（无 for 循环）
% i_idx = repmat((1:opt.num_classes)', 1, opt.num_classes); % 行索引
% j_idx = i_idx'; % 列索引
% mu_diff = mu(i_idx, :) - mu(j_idx, :); % 均值差值 (num_classes × num_classes × bands)
% sigma_avg = 0.5 * (sigma(i_idx, :) + sigma(j_idx, :)); % 平均方差
% 
% B = 0.125 * (mu_diff.^2 ./ sigma_avg) + 0.5 * log(sigma_avg ./ sqrt(sigma(i_idx, :) .* sigma(j_idx, :)));
% JM_matrix = 2 * (1 - exp(-sum(B, 3))); % 计算 JM 距离矩阵
% JMtotal = sum(JM_matrix(:)); 





%MMI
% idx=find(bX==1);
% selected_MI = opt.matrix(idx, idx);
% upper_tri = triu(selected_MI, 1);
% mmi = sum(upper_tri(:));
idx=bX==1;
selected_MI = opt.matrix(idx, idx);
mmi = sum(selected_MI,'all')/2;



%{
信噪比
[cor,~] = get_dataset_name(2);
[H, W, bands] = size(cor);
selected_bands=find(X);
nB=length(selected_bands);
cors=cor(:,:,selected_bands);    
snr_db = zeros(1, nB);
signal = zeros(1, nB);
noise = zeros(1, nB);
signal = mean(reshape(cors, [], nB), 1);      % 各波段均值
noise = std(reshape(cors, [], nB), 0, 1);    % 各波段标准差
% 噪声下限处理
noise = max(noise, eps);
% 计算SNR(dB)
snr_db = 20 * log10(signal ./ noise);
snr_db1=sum(snr_db);
%}








 %{                                     
[cor,~] = get_dataset_name(2);      波段方差
[H, W, ~] = size(cor);
selected_bands=find(X);
nB=length(selected_bands);
cors=cor(:,:,selected_bands);
band_variances = zeros(1, nB);
reshaped_data = reshape(cors, [], nB);
band_variances = var(reshaped_data, 1, 1, 'omitnan'); 
total_variance = sum(band_variances);
 %}

           





%熵

one_indices=find(bX==1);     %%f2
sum_counts=sum(opt.hist(:,one_indices),2);     
p=sum_counts./(opt.h*opt.w*60);

p=-p.*log2(p);
p(isnan(p))=0;
h=sum(p)/sum(bX);      %%f2

%所选波段数量通过a进行调节
% a=opt.select_para;                      %%f3
% s=abs((dim-a*sum(bX))/dim);          %%f3

% a=opt.band_num;                 %%f3
% s=abs(sum(bX)-a)*100;          %%f3


fitness=-sum(opt.LD(bX==1))/sum(bX)-h+mmi;   %-sum(opt.LD(bX==1))/sum(bX)
% disp(sum(bX));



end

