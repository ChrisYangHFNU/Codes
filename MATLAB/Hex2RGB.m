function c = Hex2RGB(str)
% 十六进制颜色转换成RGB颜色数组

clist = '0123456789ABCDEF';
nums = zeros(1, 6);
for i = 1:6
    nums(i) = find(str(i + 1) == clist);
end
c = zeros(1, 3);
for i = 1:3
    c(i) = 16 * (nums(2*i-1) - 1) + (nums(2*i) - 1);
end
end