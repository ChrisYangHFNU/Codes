function fx=ZAobjFun(x,k1,k2)
% ZAobjFun  ����һ��������ѧ���������������ֺ��������
% x			�Ǻ�������
% k1, k2	�Ƕ�����庯������Ĳ���
% ZZYд�� 2014/5/8

fx=k1*sin(tan(x))-k2*tan(sin(x-0.7));
