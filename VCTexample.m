warning('off')
num_vertices = 6;
angles = (0:num_vertices-1) * (2*pi / num_vertices);
x_coords = 2 * cos(angles);
y_coords = 2 * sin(angles);
A = x_coords(:);        % A是以原点为圆心的，边长为2的正六边形边界的横坐标
B = y_coords(:);        % B是以原点为圆心的，边长为2的正六边形边界的纵坐标
x=rand(10,1);           % 输入Vofs函数的十个随机点的横坐标（在A,B所围区域范围内）
y=rand(10,1);           % 输入Vofs函数的十个随机点的纵坐标（在A,B所围区域范围内）
r=10^-5;                % 收敛程度，越小越接近VCT，但0.01级别就已经精度很高，且符合直觉了
[C,D]=Vofs(A,B,x,y,r)   % Vofs函数，输出C为上次迭代的维诺质心镶嵌的质心位置，D为每个质心细胞对应的顶点位置。
voronoi(C(:,1),C(:,2)); % 根据质心绘制维诺图
axis equal      
hold on
A=[A;A(1,1)];           % 为首尾相接补充A
B=[B;B(1,1)];           % 为首尾相接补充B
plot(A,B);              % 绘制边界