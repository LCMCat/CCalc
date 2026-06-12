# CCalc

一个使用 C++17 编写的高精度命令行高级计算器，零外部依赖。

所有核心算术运算（BigInt、BigRat、BigFloat）均从零开始实现，支持 **100+ 位**有符号数和**精确符号计算**（例如 `sin(pi/2) = 1`、`ln(e) = 1`）。

## 功能特性

### 算术运算
- `+  -  *  /  ^  %` 任意精度运算
- 精确有理数运算：`10/3` → `10/3 (= 3.333...)`
- 隐式乘法：`2pi`、`3(1+2)`、`(2)pi`
- 逗号分隔表达式：`x=1, x+1, 2/2` → 分别计算并输出每个结果

### 常量
- `pi` — π，具有 1000+ 位精度
- `e` — 欧拉数，具有 1000+ 位精度

### 变量
- `x`, `y`, `A`, `B`, `C`, `D` — 用户可赋值变量（默认值：0）
- 赋值：`x=1`、`x=sin(pi/2)`、`A=10`
- 在表达式中使用：当 `x=pi` 时，`sin(x/2)` → `1`

### 自定义函数
- 定义：`f(x) := x^2+1`、`g(x,y) := x*y+1`
- 调用：`f(3)` → `10`、`g(2,5)` → `11`
- 支持所有内置函数：`h(x) := sin(x)+cos(x)`

### 向量
- `VerA`, `VerB`, `VerC`, `VerD` — 用户可赋值向量变量（默认值：空）
- 赋值：`VerA=(1,2,3)`、`VerB=(1,2)`
- 向量字面量：`(1,2,3)` — 三维向量，`(3,4)` — 二维向量
- 算术运算：`VerA+VerB`、`VerA-VerB`、`3*VerA`、`-VerA`（加减需要相同维度）
- 运算：
  - `vecmod(v)` — 模长：`vecmod((3,4))` → `5`
  - `dot(a, b)` — 点积（需要相同维度）
  - `cross(a, b)` — 叉积（仅三维）
  - `scalarmul(s, v)` — 标量乘法（也可用 `s*v` 或 `v*s`）
  - `mixed(a, b, c)` — 混合积 a·(b×c)（仅三维）
  - `proj(a, b)` — a 在 b 上的投影
  - `decompose(a, b, c)` — 二维分解：a = αb + βc，返回 (α, β)
  - `decompose(a, b, c, d)` — 三维分解：a = αb + βc + γd，返回 (α, β, γ)

### 进制转换
- 输入前缀：`0b1010`（二进制）、`0o77`（八进制）、`0xFF`（十六进制）
- 输出进制：`base 2` / `base 8` / `base 16` / `base N`（2–36）
- 直接进制运算：`0b101+0b110` → `11`
- 示例：`base 16; 255` → `FF`

### 三角函数
- `sin`, `cos`, `tan` — 特殊角度具有精确值
- `asin`, `acos`, `atan` — 反三角函数
- `sinh`, `cosh`, `tanh` — 双曲函数
- 角度模式：`deg` / `rad`

### 根与幂
- `sqrt(x)`, `cbrt(x)`, `nrt(n, x)` — 平方根/立方根/n次根
- `x^n` — 幂运算（整数指数为精确值）
- `n!` 或 `fact(n)` — 阶乘

### 对数
- `ln(x)` — 自然对数（`ln(e) = 1`）
- `lg(x)` — 常用对数（`lg(100) = 2`）
- `log(base, x)` — 任意底对数（`log(2, 8) = 3`）

### 数论
- `factor(n)` — 质因数分解（`factor(360)` → `2^3 * 3^2 * 5`）
- `euler_phi(n)` 或 `phi(n)` — 欧拉函数（`euler_phi(10)` → `4`）
- `gcd(a, b)`, `lcm(a, b)` — 最大公约数、最小公倍数
- `P(n, k)` / `perm(n, k)` — 排列数
- `C(n, k)` / `comb(n, k)` — 组合数

### 微积分
- `int(expr, var, a, b)` — 定积分（辛普森法则）
- `integrate(expr, var)` — **符号**不定积分：`integrate(x^2, x)` → `x^3/3 + C`
- `diff(expr, var)` — **符号**求导：`diff(x^2, x)` → `2*x`
- `partial(expr, var)` — **符号**偏导：`partial(x^2*y, y)` → `x^2`
- `diff(expr, var, point)` — 在某点的数值导数
- `taylor(expr, var, point, n)` — 泰勒展开：`taylor(sin(x), x, 0, 5)` → `x - x³/6 + x⁵/120`
- `limit(expr, var, val)` — 数值极限：`limit(sin(x)/x, x, 0)` → `1`
- `sum(expr, var, from, to)` — 求和
- `prod(expr, var, from, to)` — 求积
- `expand(expr)` — 多项式展开：`expand((x+1)^3)` → `1+3*x+3*x^2+x^3`

### 积分表
- `inttable(expr)` — 查询常见函数的原函数
- 支持：x^n、sin(x)、cos(x)、exp(x)、ln(x)、1/x 及其通过和差法则的组合

### 复数
- `complex(a, b)` — 构造 `a + bi`
- `re(z)`, `im(z)`, `conj(z)`, `arg(z)`, `mod(z)` — 实部、虚部、共轭、幅角、模
- 自动复数化：`sqrt(-1)` → `i`、`ln(-1)` → `pi*i`

### 其他
- `abs(x)`, `floor(x)`, `ceil(x)`, `round(x)`, `sign(x)` — 绝对值、向下取整、向上取整、四舍五入、符号
- `max(a, b)`, `min(a, b)` — 最大值、最小值
- `rand()` — [0, 1) 范围随机数
- `randint(a, b)` — 随机整数
- `convert(val, from, to)` — 单位转换（长度、质量、时间、面积、体积、速度、数据、温度）

### 概率
- `normal_pdf(x, mu, sigma)` — 正态分布概率密度
- `normal_cdf(x, mu, sigma)` — 正态分布累积分布

### 几何
- `triangle_area(a, b, c)` — 根据边长计算三角形面积（海伦公式）：`triangle_area(3, 4, 5)` → `6`
- `polygon_area(points)` — 根据顶点计算多边形面积（鞋带公式）：`polygon_area([0,0;1,0;1,1;0,1])` → `1`
- `distance(p1, p2)` — 两点间欧几里得距离：`distance([0,0], [3,4])` → `5`
- `line_intersect(y=expr1, y=expr2)` — 两直线交点：`line_intersect(y=x, y=2-x)` → `(1, 1)`

### 方程求解
- `solve(a, b, c)` — 二次方程 ax²+bx+c=0（**精确根式解**）
- `solve(a, b, c, d)` — 三次方程 ax³+bx²+cx+d=0（卡尔达诺公式）
- `solve(a, b, c, d, e)` — 四次方程 ax⁴+bx³+cx²+dx+e=0（费拉里方法）
- 支持复根：`solve(1,0,1)` → `(i, -i)`
- 包含判别式、局部/全局极值及极值点

### 矩阵运算
- 矩阵字面量：`[1,2;3,4]`（行用 `;` 分隔，列用 `,` 分隔）
- 算术运算：`+`、`-`、`*`（矩阵乘法）、标量乘法
- `det(m)` — 行列式（精确有理数运算）
- `inv(m)` — 逆矩阵（精确有理数运算）
- `eigen(m)` — 特征值（2×2 和 3×3）
- `trace(m)` — 迹
- `transpose(m)` — 转置
- `identity(n)` 或 `eye(n)` — n×n 单位矩阵
- **美观输出**，使用 Unicode 制表符

### 统计
- `mean(a,b,c,...)` — 算术平均值
- `stddev(a,b,c,...)` — 样本标准差
- `variance(a,b,c,...)` — 样本方差
- `median(a,b,c,...)` — 中位数

### 表达式化简
- `simplify(expr)` — 计算并显示化简结果
- 示例：`simplify(6/3)` → `6/3 = 2`

### 递推序列
- `recur(expr, var, init, n)` — 计算递推数列的第 n 项
- 示例：`recur(a+1, a, 1, 5)` → `6`（等差数列：1,2,3,4,5,6）
- 示例：`recur(a*2, a, 1, 10)` → `1024`（等比数列：1,2,4,...,1024）

### 函数值表
- `table(expr, var, from, to[, step])` — 输出函数值表
- 示例：`table(x^2, x, 0, 5)` → x = 0,1,2,3,4,5 时 x² 的值
- 示例：`table(sin(x), x, 0, 6, 1)` → x = 0,1,...,6 时 sin(x) 的值

### 拉格朗日乘数法
- `lagrange(f, g, x, y)` — 在约束 g(x,y)=0 下求 f(x,y) 的极值
- 输出偏导数、方程组和临界点
- 示例：`lagrange(x+y, x^2+y^2-1, x, y)` → 在单位圆上求 x+y 的最大/最小值

### LaTeX 输出
- `latex` 命令切换 LaTeX 格式化模式
- 分数：`1/2` → `\frac{1}{2}`
- 平方根：`sqrt(2)` → `\sqrt{2}`
- 幂：`x^2` → `x^{2}`
- 矩阵：`\begin{pmatrix}1 & 2 \\ 3 & 4\end{pmatrix}`
- 希腊字母：`pi` → `\pi`

### 批处理模式
- `ccalc -f file.txt` — 从文件读取表达式，逐行输出结果
- 以 `#` 开头的行是注释
- 变量赋值在批处理模式下有效

### REPL 功能
- **历史记录导航**：上/下箭头键浏览历史输入
- **Tab 补全**：输入 `si` + Tab → `sin(`
- **语法高亮**：数字（灰色）、函数（绿色）、运算符（黄色）
- **多行输入**：行末加 `\` 继续，或括号未闭合时自动继续
- **颜色主题**：`theme` 列出主题，`theme dark`/`theme neon`/`theme light`/`theme none` 切换主题
- **中文模式**：`zh` 切换到中文，`en` 切换回英文
- **配置文件**：设置（精度、角度模式、主题、自定义函数）退出时自动保存

### 函数绘图
- `graph(expr)` — 绘制 y=expr（例如 `graph(sin(x))`、`graph(x^2-3*x+1)`）
- `graph_implicit(expr)` — 绘制 f(x,y)=0（例如 `graph_implicit(x^2+y^2-1)` 绘制圆）
- `graph_param(t, x_expr, y_expr)` — 参数曲线（例如 `graph_param(t, cos(t), sin(t))` 绘制圆）
- `graph_polar(r_expr)` — 极坐标曲线 r=f(θ)（例如 `graph_polar(1+cos(theta))` 绘制心形线）
- `graph3d(expr)` — 三维曲面图 z=f(x,y)（例如 `graph3d(sin(x)*cos(y))`）
- 基于 Dear ImGui + DirectX 11
- 或直接运行：`ccalc_graph "sin(x)"`、`ccalc_graph -i "x^2+y^2-1"`、`ccalc_graph -p "cos(t)" "sin(t)"`、`ccalc_graph -l "1+cos(theta)"`
- 或直接运行三维：`ccalc_graph3d "sin(x)*cos(y)"`
- 简洁界面：默认只显示图形
  - 点击 **Settings >>** 按钮切换设置面板
  - 设置：添加/删除函数、图形类型、视图范围、网格、角度模式、参数范围
- 功能：
  - 同时绘制最多 8 个函数，不同颜色
  - 四种图形类型：显式 y=f(x)、隐式 f(x,y)=0、参数式、极坐标
  - 隐式使用 Marching Squares 算法绘制平滑等高线
  - 缩放（滚轮）、平移（左键拖拽）、悬停显示坐标
  - 1:1 纵横比保证形状准确
  - 自动 Y 范围、网格、坐标轴显示
  - 角度模式（弧度/角度）
  - 表达式语法与计算器相同

## 精确值示例

| 输入 | 输出 |
|------|------|
| `sin(pi/2)` | `1` |
| `cos(pi/3)` | `1/2 (= 0.5)` |
| `tan(pi/4)` | `1` |
| `sin(pi/6)` | `1/2 (= 0.5)` |
| `cos(pi/4)` | `sqrt(2)/2 ~= 0.7071...` |
| `sin(pi/3)` | `sqrt(3)/2 ~= 0.8660...` |
| `asin(1)` | `pi/2 ~= 1.5707...` |
| `acos(1/2)` | `pi/3 ~= 1.0471...` |
| `atan(1)` | `pi/4 ~= 0.7853...` |
| `ln(e)` | `1` |
| `lg(100)` | `2` |
| `log(2, 8)` | `3` |
| `sqrt(4)` | `2` |
| `cbrt(27)` | `3` |
| `factor(360)` | `2^3 * 3^2 * 5` |
| `sum(k^2, k, 1, 10)` | `385` |
| `deg; sin(90)` | `1` |
| `deg; cos(180)` | `-1` |
| `x=pi; sin(x/2)` | `1` |
| `0xFF` | `255` |
| `0b1010` | `10` |
| `base 16; 255` | `FF` |
| `VerA=(1,2,3)` | `(1, 2, 3)` |
| `vecmod((3,4))` | `5` |
| `dot((1,2),(3,4))` | `11` |
| `cross((1,0,0),(0,1,0))` | `(0, 0, 1)` |
| `decompose((1,1),(1,0),(0,1))` | `(1, 1)` |
| `x=1, x+1, 2/2, 1-2*5` | `1`, `2`, `1`, `-9` |
| `solve(1,2,1)` | `(-1, -1)` |
| `solve(1,-1,-1)` | `(sqrt(5)/2+1/2, -sqrt(5)/2+1/2)` |
| `solve(1,0,0,-1)` | `(1, -0.5+0.866i, -0.5-0.866i)` |
| `solve(1,0,1)` | `(i, -i)` |
| `det([1,2;3,4])` | `-2` |
| `inv([1,2;3,4])` | `[[-2, 1]; [3/2, -1/2]]` |
| `[1,2;3,4]*[5;6]` | `[[17]; [39]]` |
| `mean(1,2,3,4,5)` | `3` |
| `stddev(1,2,3,4,5)` | `sqrt(2) ~= 1.414...` |
| `median(1,2,3,4,5)` | `3` |
| `simplify(6/3)` | `6/3 = 2` |
| `diff(x^2, x)` | `2*x` |
| `integrate(x^2, x)` | `x^3/3 + C` |
| `partial(x^2*y, y)` | `x^2` |
| `expand((x+1)^3)` | `1+3*x+3*x^2+x^3` |
| `normal_pdf(0, 0, 1)` | `0.3989...` |
| `triangle_area(3, 4, 5)` | `6` |
| `distance([0,0], [3,4])` | `5` |
| `taylor(sin(x), x, 0, 5)` | `x - x³/6 + x⁵/120` |
| `limit(sin(x)/x, x, 0)` | `1` |
| `f(x) := x^2+1; f(3)` | `10` |
| `recur(a+1, a, 1, 5)` | `6` |
| `euler_phi(10)` | `4` |
| `phi(30)` | `8` |

## 构建

### 前置要求
- C++17 兼容编译器（GCC ≥ 9、Clang ≥ 10、MSVC ≥ 2019）
- Make
- DirectX 11（用于绘图窗口，仅 Windows）

### 编译

```bash
make            # 构建 ccalc.exe、ccalc_graph.exe 和 ccalc_graph3d.exe
make ccalc.exe  # 仅构建计算器
make ccalc_graph.exe  # 仅构建 2D 绘图窗口
make ccalc_graph3d.exe  # 仅构建 3D 绘图窗口
```

### 清理

```bash
make clean
```

## 使用方法

运行计算器：

```bash
./ccalc
```

### 批处理模式

```bash
./ccalc -f expressions.txt
```

### 命令

| 命令 | 描述 |
|------|------|
| `deg` | 切换到角度模式 |
| `rad` | 切换到弧度模式 |
| `base N` | 设置输出进制（2–36） |
| `prec N` | 设置精度为 N 位（1–10000） |
| `latex` | 切换 LaTeX 输出模式 |
| `zh` | 切换到中文模式 |
| `en` | 切换到英文模式 |
| `theme [name]` | 列出或切换颜色主题 |
| `ans` | 上一个答案 |
| `graph(expr)` | 绘制函数 y=expr |
| `graph3d(expr)` | 三维曲面图 z=expr |
| `help` | 显示帮助 |
| `quit` / `exit` | 退出 |

### 交互会话

```
CCalc v2.0 - Command-line Advanced Calculator
Type 'help' for help, 'quit' to exit.
Up/Down: history, Tab: complete, latex: toggle LaTeX mode

CCalc> diff(x^3, x)
3*x^2
CCalc> taylor(exp(x), x, 0, 4)
1 + 1*x + 1/2*x^2 + 1/6*x^3 + 1/24*x^4
CCalc> f(x) := x^2+1
Defined: f(x) = x^2+1
CCalc> f(3)
10
CCalc> limit(sin(x)/x, x, 0)
1
CCalc> latex
LaTeX mode: ON
CCalc> 1/2
\frac{1}{2}
CCalc> det([1,2;3,4])
-2
CCalc> quit
Goodbye!
```

## 架构

```
┌──────────┐    ┌─────────┐    ┌────────┐    ┌──────────┐    ┌───────┐
│  输入    │───>│  词法器  │───>│ 解析器 │───>│  求值器  │───>│ 值    │───> 输出
│ (字符串)  │    │ (词元)   │    │  (AST) │    │  (求值)  │    │(结果) │
└──────────┘    └─────────┘    └────────┘    └──────────┘    └───────┘
```

### 核心类型

| 类型 | 描述 |
|------|------|
| `BigInt` | 任意精度整数（基数为 10^9） |
| `BigRat` | 任意精度有理数（自动约分） |
| `BigFloat` | 任意精度浮点数（科学计数法） |
| `SurdsExpr` | 精确符号表达式（根式 + π + e） |
| `Value` | 判别联合：SURDS / FLOAT / COMPLEX / VECTOR / MATRIX / STRING / ERROR |

### 关键设计决策

- **SurdsExpr** 支持精确符号计算。特殊根式值：`-1` = π、`-2` = e、`≥2` = √n。
- **精确优先求值**：运算尽可能保持在精确（SurdsExpr）域，仅在必要时回退到 BigFloat。
- **无外部依赖**：所有高精度算术运算均从零实现。
- **硬编码常量**：π、e、ln(2) 存储 1000+ 位以保证最大精度。
- **符号微分**：AST 变换方法，支持链式法则、乘积法则、商法则。

## 许可证

MIT