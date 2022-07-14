function c = RGB2MatlabColor(rgb)
% RGB转MATLAB颜色向量

c = round(100 * rgb/255) / 100;
end