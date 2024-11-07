

# Dot line 


1. 以x,y为变量，进行线性回归

2. 获取回归方差的系数 （β, b）

3. 计算所有值到直线的垂直距离： |A\*x + B\*y + C| \/ sqrt(A^2 + B ^2)，特殊形式： |kx+b-y|\/sqrt(1+k^2)

4. 获取所有距离的SD，以 N\*SD 在回归曲线上下进行可视化和点的筛选

<image src="./comparisions_dotplot.jpg">