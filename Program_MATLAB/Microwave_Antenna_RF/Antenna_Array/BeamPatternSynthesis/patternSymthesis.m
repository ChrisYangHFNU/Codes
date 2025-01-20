clear 

options = gaoptimset('PopulationSize',30,... %种群包含个体数目
    'EliteCount',15,...  % 种群中精英个体数目
    'CrossoverFraction',0.8,... % 交叉后代比例
    'MigrationFraction',0.2,... % 变异概率
    'Generations',2000,... % 迭代代数
    'StallGenLimit',2000,... % 停止代数
    'TolFun',1e-100,... % 适应度函数偏差
    'PlotFcns',{@gaplotbestf,@gaplotbestindiv,@gaplotstopping}); % 绘制最优个体适应度函数与最优个体

% options = gaoptimset('ParetoFraction',0.3,...
%     'PopulationSize',40,...
%     'Generations',200,...
%     'StallGenLimit',200,...
%     'TolFun',1e-100,...
%     'PlotFcns',@gaplotpareto);

eleNum = 20;
phaLB = -180*ones(1,eleNum);ampLB = zeros(1,eleNum);
phaUB = 180*ones(1,eleNum);ampUB = ones(1,eleNum);
LB = [phaLB ampLB];
UB = [phaUB ampUB];
% 运行遗传算法
[x,fval,exitFlag,Output] = ga(@objFun_beamshape,2*eleNum,[],[],[],[],LB,UB,[],[],options)
% [x,fval,exitFlag,output] = gamultiobj(@objFun_beamshape,2*eleNum,[],[],[],[],LB,UB,options);
