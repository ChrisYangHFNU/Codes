function c = RGB2MatlabColor(rgb)
% RGBתMATLAB��ɫ����

c = round(100 * rgb/255) / 100;
end