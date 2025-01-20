function fx=ZAobjFun(x,k1,k2)
% ZAobjFun  表述一个含参数学函数，供创建各种函数句柄用
% x			是函数变量
% k1, k2	是定义具体函数所需的参数
% ZZY写于 2014/5/8

fx=k1*sin(tan(x))-k2*tan(sin(x-0.7));
