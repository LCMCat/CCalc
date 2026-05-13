# CCalc 技术详细报告

## 1. 项目概述

CCalc 是一个用 C++17 实现的命令行高级计算器，支持高精度运算（100位有符号数）、精确值计算（如 `sin(pi/2) = 1`）以及丰富的数学函数库。项目不依赖任何外部高精度数学库（如 GMP），所有核心运算从零实现。

## 2. 系统架构

### 2.1 模块结构

```
CCalc/
├── ccalc.h          # 统一头文件：所有类声明、接口定义
├── bignum.cpp       # 高精度数值运算：BigInt, BigRat, BigFloat
├── value.cpp        # 值类型系统：SurdsExpr, Value, ComplexVal
├── parser.cpp       # 词法分析与语法分析：Lexer, Parser
├── evaluator.cpp    # 表达式求值与函数实现：Evaluator
├── main.cpp         # REPL 交互界面
├── graph.cpp        # 函数图像绘制窗口（ImGui + DirectX 11）
├── Makefile         # 构建脚本
└── imgui-1.92.7/    # Dear ImGui 库（仅 graph 使用）
```

### 2.2 数据流

```
用户输入 → Lexer(词法分析) → Token流 → Parser(语法分析) → AST → Evaluator(求值) → Value → 格式化输出
```

## 3. 核心类设计

### 3.1 BigInt — 任意精度整数

**存储方式**：采用 `BASE = 10^9` 的进位制，`vector<int64_t>` 存储各位数字（小端序），`bool neg_` 标记符号。

```
数值 = (-1)^neg_ * Σ digits_[i] * BASE^i
```

**关键算法**：

| 运算 | 算法 | 时间复杂度 |
|------|------|-----------|
| 加/减 | 逐位运算 + 进位传播 | O(n) |
| 乘法 | Karatsuba（递归分治），n<64 回退学校乘法 | O(n^1.585) |
| 除法 | 长除法 | O(n²) |
| GCD | 辗转相除法 | O(n²) |
| 幂运算 | 快速幂 | O(log n) 次乘法 |
| 阶乘 | 逐次乘法 | O(n²) |

**特殊设计**：
- `divmod()` 同时返回商和余数，避免重复计算
- Karatsuba 乘法：当操作数位数 ≥ 64 时自动切换到 Karatsuba 分治算法
  - 将 x = x₁·B^m + x₀, y = y₁·B^m + y₀
  - 递归计算 z₀ = x₀y₀, z₂ = x₁y₁, z₁ = (x₀+x₁)(y₀+y₁) - z₀ - z₂
  - 结果 = z₂·B^(2m) + z₁·B^m + z₀
  - 辅助函数：`add_vec()`（向量加法）、`sub_vec()`（向量减法）、`shift()`（左移偏移）
  - 递归基（n < 32）使用 O(n²) 学校乘法
- `to_int64()` 提供与原生类型的互转
- `pow()` 使用二进制快速幂算法
- `to_base_string(base)` 将整数转换为 2~36 进制字符串
- `from_base_string(s, base)` 从 2~36 进制字符串构造整数

### 3.2 BigRat — 任意精度有理数

**存储方式**：`BigInt num_`（分子）+ `BigInt den_`（分母），始终保证 `den_ > 0` 且最简（GCD = 1）。

**自动简化**：每次构造和运算后调用 `simplify()`，通过 `BigInt::gcd()` 约分。

**输出格式**：
- `to_string()`：分数形式 `a/b` 或整数形式 `a`
- `to_decimal_string(prec)`：小数展开，精确到指定位数

### 3.3 BigFloat — 任意精度浮点数

**存储方式**：科学记数法 `mantissa_ × 10^exp_`，其中 `mantissa_` 是 BigInt。

```
数值 = mantissa_ * 10^(exp_ - mantissa_digits + 1)
```

**归一化**：`normalize(prec)` 去除尾随零，并按精度截断多余位数。

**关键数学函数实现**：

| 函数 | 算法 | 说明 |
|------|------|------|
| `pi()` | 硬编码 1000+ 位 | 首次调用后缓存 |
| `e_val()` | 硬编码 1000+ 位 | — |
| `ln2()` | 硬编码 1000+ 位 | — |
| `exp(x)` | Taylor 级数 Σ x^n/n! | 收敛快 |
| `ln(x)` | AGM 算法 | 比级数收敛更快 |
| `sin/cos` | Taylor 级数 + 参数归约 | 先归约到 [0, π/4] |
| `asin/acos` | Taylor + 恒等变换 | |x|>0.5 时用互补公式 |
| `atan` | Taylor 级数 | |x|>1 时用倒数公式 |
| `sinh/cosh/tanh` | 定义式 | 基于 exp |
| `sqrt` | Newton 迭代 | 二次收敛 |
| `cbrt` | Newton 迭代 | 三次收敛 |
| `nrt` | Newton 迭代 | n 次收敛 |
| `pow` | exp(y*ln(x)) | 整数幂优化为连乘 |

**字符串构造函数**支持：
- 普通小数：`"3.14159"`
- 科学记数法：`"1e-120"`, `"2.5E100"`

### 3.4 SurdsExpr — 根式精确表达式

**设计动机**：许多数学结果无法用有理数精确表示，但可以用根式的线性组合表示。例如 `cos(pi/4) = sqrt(2)/2`。

**存储方式**：`vector<SurdsTerm> terms`，每个 `SurdsTerm` 包含：
- `BigRat coeff`：系数
- `BigInt radicand`：被开方数

**特殊 radicand 约定**：

| radicand | 含义 | 示例 |
|----------|------|------|
| 1 | 有理数项 | `3/4` → `SurdsTerm(3/4, 1)` |
| -1 | π 常数 | `pi/2` → `SurdsTerm(1/2, -1)` |
| -2 | e 常数 | `3*e` → `SurdsTerm(3, -2)` |
| ≥2 | √radicand | `2*sqrt(3)` → `SurdsTerm(2, 3)` |

**运算规则**：
- **加法**：合并同类项（相同 radicand）
- **乘法**：交叉相乘，`sqrt(a)*sqrt(b) = sqrt(a*b)`，自动做无平方因子分解
- **π², e², π*e** 等超越数乘积回退到浮点数
- `make_square_free()` 将 radicand 分解为 `k² * r`，提取 `k` 到系数中

**显示格式**：
- `sqrt(2)/2` 而非 `1/2*sqrt(2)`
- `3*pi/4` 而非 `3/4*pi`
- `2*e` 而非 `2*sqrt(-2)`

### 3.5 Value — 统一值类型

**类型枚举**：

```cpp
enum Type { SURDS, FLOAT, COMPLEX, VECTOR, MATRIX, STRING, ERROR };
```

| 类型 | 用途 | 存储字段 |
|------|------|---------|
| SURDS | 精确值（有理数/根式/π/e） | `surds: SurdsExpr` |
| FLOAT | 浮点近似值 | `float_val: BigFloat` |
| COMPLEX | 复数 | `complex: ComplexVal` |
| VECTOR | 向量 | `vec: vector<Value>` |
| MATRIX | 矩阵 | `mat: vector<vector<Value>>` |
| STRING | 字符串（因数分解/求导/积分表等结果） | `error_msg: string` |
| ERROR | 错误信息 | `error_msg: string` |

**MATRIX 类型**：

- `mat` 字段存储 `vector<vector<Value>>`，每行一个 `vector<Value>`
- `mat_rows()` / `mat_cols()` 返回矩阵维度
- `make_matrix(rows)` 静态工厂方法构造矩阵
- `is_matrix()` 判断是否为矩阵类型

**类型提升规则**：
- SURDS ⊕ FLOAT → FLOAT
- SURDS/FLOAT ⊕ COMPLEX → COMPLEX
- VECTOR ⊕ VECTOR → VECTOR（同维度加减、标量乘法）
- MATRIX ⊕ MATRIX → MATRIX（同维度加减、矩阵乘法）
- MATRIX ⊕ 标量 → MATRIX（标量乘法）
- 任何 ⊕ ERROR → ERROR

**ComplexVal**：`shared_ptr<Value> real, imag`，支持任意精度的实部和虚部。

## 4. 解析器设计

### 4.1 Lexer — 词法分析器

**Token 类型**：

```
NUMBER, IDENTIFIER, PLUS, MINUS, STAR, SLASH, CARET, PERCENT,
LPAREN, RPAREN, LBRACKET, RBRACKET, COMMA, SEMICOLON, BANG, EQUAL, LT, GT, LE, GE, NEQ,
COLON_EQUAL,
END_OF_INPUT, ERROR
```

**新增 Token 说明**：

| Token | 字符 | 用途 |
|-------|------|------|
| LBRACKET | `[` | 矩阵字面量左括号 |
| RBRACKET | `]` | 矩阵字面量右括号 |
| SEMICOLON | `;` | 矩阵行分隔符 |
| COLON_EQUAL | `:=` | 自定义函数定义 |

**隐式乘法**：自动在以下相邻 Token 间插入 `*`：
- 数字后接标识符：`2pi` → `2*pi`
- 数字后接左括号：`3(1+2)` → `3*(1+2)`
- 右括号后接数字：`(2)3` → `(2)*3`
- 右括号后接左括号：`(1)(2)` → `(1)*(2)`
- 右括号后接标识符：`(2)pi` → `(2)*pi`

**进制前缀**：Lexer 支持以下前缀将非十进制输入转为十进制数值：
- `0b` / `0B`：二进制，如 `0b1010` → 10
- `0o` / `0O`：八进制，如 `0o77` → 63
- `0x` / `0X`：十六进制，如 `0xFF` → 255

前缀后的字符必须属于对应进制的合法数字集，否则解析在该字符处停止。

### 4.2 Parser — 语法分析器

**递归下降解析**，运算符优先级（由低到高）：

| 优先级 | 运算符 | 结合性 |
|--------|--------|--------|
| 1 | `< > <= >=` | 左结合 |
| 2 | `+ -` | 左结合 |
| 3 | `* / %` | 左结合 |
| 4 | `^` | 右结合 |
| 5 | 一元 `-` | 前缀 |
| 6 | `!` (阶乘) | 后缀 |

**AST 节点类型**：

```
NUMBER      — 数值字面量 (BigRat)
CONSTANT    — 常量 (pi, e)
VARIABLE    — 变量 (x, k, ans, VerA, VerB, VerC, VerD)
BINOP       — 二元运算 (+, -, *, /, ^, %)
UNARYOP     — 一元运算 (-)
FUNCTION    — 函数调用 (sin, cos, ...)
FACTORIAL   — 阶乘后缀 (n!)
VEC_LITERAL — 向量字面量 ((1,2,3))
MAT_LITERAL — 矩阵字面量 ([1,2;3,4])
```

**矩阵字面量解析**：Parser 的 `primary()` 函数检测 `[` 开头的语法，解析为 `MAT_LITERAL` AST 节点：

- `[` 开始矩阵字面量
- 逗号 `,` 分隔同行元素
- 分号 `;` 分隔行
- `]` 结束矩阵字面量
- 示例：`[1,2;3,4]` → 2×2 矩阵，`[1,2,3]` → 1×3 矩阵（行向量）
- `mat_rows` 字段存储 `vector<vector<ASTPtr>>`

**已知函数列表**：60+ 个内置函数名，用于区分函数调用和变量名乘法。

## 5. 求值器设计

### 5.1 精确值计算策略

**核心原则**：尽可能在 SurdsExpr 域中运算，只在无法精确表示时回退到 BigFloat。

**三角函数精确值**：`try_exact_trig(pi_coeff, func)` 查表返回 24 个特殊角度的精确值：

| π 的系数 | sin | cos | tan |
|----------|-----|-----|-----|
| 0 | 0 | 1 | 0 |
| 1/6 | 1/2 | √3/2 | √3/3 |
| 1/4 | √2/2 | √2/2 | 1 |
| 1/3 | √3/2 | 1/2 | √3 |
| 1/2 | 1 | 0 | ∞ |
| 2/3 | √3/2 | -1/2 | -√3 |
| 3/4 | √2/2 | -√2/2 | -1 |
| 5/6 | 1/2 | -√3/2 | -√3/3 |
| 1 | 0 | -1 | 0 |
| 7/6 ~ 11/6 | 对称值 | | |

**角度模式处理**：
- 弧度模式：直接检查输入是否为 π 的有理数倍
- 度数模式：`to_radians()` 将度数转为 `SurdsExpr(deg/180, -1)` 即 π 系数形式，然后走同一套精确值检测

**对数精确值**：
- `ln(e) = 1`, `ln(1) = 0`, `ln(-1) = πi`
- `lg(10^n) = n`, `lg(100) = 2`
- `log(a, a^n) = n`, `log(e, e) = 1`

**幂运算精确值**：
- 有理数底有理数整数幂 → `BigRat::pow()`
- BigFloat 底整数幂 → 连乘（避免 `exp(n*ln(x))` 的精度损失）
- SurdsExpr 底整数幂 → 连乘（保持根式精确性）
- e 的整数幂 → `BigFloat::exp(n)`

### 5.2 微积分运算

**定积分** — Simpson 法则：
```
∫[a,b] f(x)dx ≈ h/3 * [f(a) + f(b) + 4*Σf(x_{2i-1}) + 2*Σf(x_{2i})]
```
- 1000 个区间
- 通过 `substitute()` 替换变量值

**数值微分** — 五点中心差分公式（O(h⁴) 精度）：
```
f'(x) ≈ [-f(x-2h) + 8f(x-h) - 8f(x+h) + f(x+2h)] / (12h)
```
- `h = 10^(-precision/3 - 5)`

**符号求导** — `diff(expr, var)`：

通过 AST 变换实现符号求导，`sdiff(node, var)` 递归遍历 AST 并返回导数的 AST：

| 节点类型 | 求导规则 | 示例 |
|----------|---------|------|
| NUMBER | d/dx(c) = 0 | `sdiff(3, x)` → `0` |
| CONSTANT | d/dx(π) = 0 | `sdiff(pi, x)` → `0` |
| VARIABLE | d/dx(x) = 1, d/dx(y) = 0 | `sdiff(x, x)` → `1` |
| BINOP `+` | (f+g)' = f'+g' | `sdiff(x+1, x)` → `1+0` |
| BINOP `-` | (f-g)' = f'-g' | |
| BINOP `*` | (fg)' = f'g + fg'（乘积法则） | `sdiff(x*sin(x), x)` → `1*sin(x)+x*cos(x)` |
| BINOP `/` | (f/g)' = (f'g - fg')/g²（商法则） | |
| BINOP `^` | 三种情况（见下表） | |
| UNARYOP `-` | (-f)' = -f' | |
| FUNCTION sin | (sin u)' = u'·cos(u) | |
| FUNCTION cos | (cos u)' = -u'·sin(u) | |
| FUNCTION tan | (tan u)' = u'·sec²(u) | |
| FUNCTION ln | (ln u)' = u'/u | |
| FUNCTION exp | (exp u)' = u'·exp(u) | |
| FUNCTION sqrt | (√u)' = u'/(2√u) | |

**幂运算求导的三种情况**（`f^g` 对 `x` 求导）：

| 情况 | 条件 | 公式 |
|------|------|------|
| 常数幂 | f 含 x，g 不含 x | f'·g·f^(g-1) |
| 常数底 | f 不含 x，g 含 x | f^g · g'·ln(f) |
| 一般幂 | f 和 g 都含 x | f^g · (g'·ln(f) + g·f'/f)（对数微分法） |

**实现细节**：
- 通过 `ast_to_string(node)` 将导数 AST 转为可读字符串输出
- `sdiff` 不做化简，保留完整的符号形式
- `diff(expr, var, point)` 三参数形式仍使用数值微分（五点差分）

**泰勒展开** — `taylor(expr, var, point, n)`：

计算 f(x) 在 x=a 处的 n 阶泰勒展开：

```
f(x) ≈ Σ_{k=0}^{n} f^(k)(a)/k! · (x-a)^k
```

**实现方式**：
1. `current_deriv` 初始为原始表达式 AST
2. 每次迭代：
   - 用 `substitute(current_deriv, var, point)` 求值 f^(k)(a)
   - 系数 = f^(k)(a) / k!
   - 用 `sdiff(current_deriv, var)` 计算下一阶导数
   - `factorial` 逐次乘以 (k+1) 维护 k!
3. 有理系数保持精确（BigRat 运算），非有理系数回退到 BigFloat
4. 输出格式：`c₀ + c₁*x + c₂*x^2 + ... + cₙ*x^n`

**极限计算** — `limit(expr, var, val)`：

数值逼近法，在目标点两侧逐渐缩小距离：

1. ε 序列：`{1e-4, 1e-6, 1e-8, 1e-10, 1e-12}`
2. 对每个 ε，计算 `f(val+ε)` 和 `f(val-ε)` 的平均值
3. 如果相邻两次结果之差 < 1e-10，返回当前值
4. 否则继续缩小 ε
5. 若所有 ε 用完仍未收敛，返回最后一次计算值

**积分表** — `inttable(expr)`：

基于 AST 模式匹配查找常见函数的不定积分：

| 表达式模式 | 不定积分 |
|-----------|---------|
| x | x²/2 |
| x^n (n≠-1) | x^(n+1)/(n+1) |
| x^(-1) | ln\|x\| |
| sin(x) | -cos(x) |
| cos(x) | sin(x) |
| exp(x) | exp(x) |
| ln(x) | x·ln(x) - x |
| 1/x | ln\|x\| |
| f(x) + g(x) | ∫f + ∫g（逐项查表） |
| f(x) - g(x) | ∫f - ∫g（逐项查表） |
| -f(x) | -(∫f) |

**实现方式**：
- 递归遍历 AST，对每个子表达式尝试匹配已知模式
- `find_var` lambda 自动检测表达式中的变量名
- `is_const` lambda 判断子表达式是否不含变量
- 无法匹配时返回 "No known antiderivative for ..."

**递推数列** — `recur(expr, var, init, n)`：

计算递推关系 a_{n} = f(a_{n-1}) 的第 n 项：

1. `current = init_val`
2. 循环 i = 1 到 n：
   - `current = substitute(expr, var, current)`
3. 返回 `current`

**示例**：
- `recur(a+1, a, 1, 5)` → 6（等差：1,2,3,4,5,6）
- `recur(a*2, a, 1, 10)` → 1024（等比：1,2,4,...,1024）
- `recur(a^2, a, 2, 3)` → 256（a₁=2, a₂=4, a₃=16, a₄=256）

### 5.3 复数运算

**表示**：`ComplexVal { real: shared_ptr<Value>, imag: shared_ptr<Value> }`

**运算**：
- 加法：`(a+bi) + (c+di) = (a+c) + (b+d)i`
- 乘法：`(a+bi)(c+di) = (ac-bd) + (ad+bc)i`
- `re(z)`, `im(z)`, `conj(z)`, `arg(z)`, `mod(z)`

**自动复数化**：`sqrt(-1) = i`, `ln(-1) = πi`

### 5.4 单位换算

**实现方式**：基于字符串的转换因子表，避免 double 精度损失。

**支持类别**：
- 长度：m, km, cm, mm, μm, nm, in, ft, yd, mi, nmi, li, au
- 质量：kg, g, mg, t, lb, oz, jin, liang
- 时间：s, ms, min, h, d, wk, yr
- 面积：m², km², cm², ha, acre, mu
- 体积：m³, L, mL, gal, qt, cup
- 速度：m/s, km/h, mph, kn, c
- 数据：B, KB, MB, GB, TB, bit
- 温度：C, F, K（大小写不敏感，非线性换算）

### 5.5 变量系统

**预定义变量**：`x`, `y`, `A`, `B`, `C`, `D`，默认值为 0。

**赋值语法**：在 REPL 中通过 `变量名=表达式` 赋值，如：
- `x=1` — 将 x 赋值为 1
- `x=sin(pi/2)` — 将 x 赋值为 sin(pi/2) 的结果（1）
- `x=pi` — 将 x 赋值为 π

**实现方式**：
1. REPL 层面检测 `=` 号，将输入分割为变量名和表达式
2. 求值表达式部分，将结果存入 `Evaluator::variables_` 映射表
3. 后续表达式中遇到该变量名时，从映射表中读取值

**变量与函数名的冲突处理**：
- `C` 和 `P` 同时是变量名和函数名（`C(n,k)` 组合，`P(n,k)` 排列）
- 赋值语法 `C=3` 优先解释为变量赋值
- 表达式中 `C(10,3)` 优先解释为函数调用，单独 `C` 解释为变量

### 5.6 进制系统

**输入进制**：Lexer 支持以下前缀将非十进制输入转为十进制数值：
- `0b` / `0B`：二进制，如 `0b1010` → 10
- `0o` / `0O`：八进制，如 `0o77` → 63
- `0x` / `0X`：十六进制，如 `0xFF` → 255

前缀后的字符必须属于对应进制的合法数字集，否则解析在该字符处停止。

**直接进制运算**：不同进制的输入可以在同一表达式中混合使用，运算在十进制下进行：
- `0b101 + 0b110` → `11`（5 + 6 = 11）
- `0xFF - 0o10` → `247`（255 - 8 = 247）
- `0b100 * 0xA` → `40`（4 * 10 = 40）

**输出进制**：通过 `base N` 命令设置（2~36），默认为 10。

**实现方式**：
- `BigInt::to_base_string(base)` 将整数转为指定进制字符串
- `BigInt::from_base_string(s, base)` 从指定进制字符串构造整数
- 字母 A~Z 表示 10~35
- 仅对整数结果应用进制转换，分数和根式仍以十进制显示

**示例**：
| 命令 | 输出 |
|------|------|
| `base 16; 255` | `FF` |
| `base 2; 10` | `1010` |
| `base 8; 64` | `100` |
| `base 4; 15` | `33` |
| `base 36; 35` | `Z` |

### 5.7 向量系统

**向量变量**：`VerA`, `VerB`, `VerC`, `VerD`，默认为空向量（0维）。

**赋值语法**：`VerA=(1,2,3)`，`VerB=(1,2)`

**存储方式**：`Value::vec` 为 `std::vector<Value>`，每个元素是一个 Value（可以是 SURDS、FLOAT 等任意类型）。

**向量字面量解析**：Parser 的 `primary()` 函数检测 `(expr, expr, ...)` 语法，当括号内第一个表达式后紧跟逗号时，解析为 `VEC_LITERAL` AST 节点，而非分组括号。

**维度约束**：
- 向量加法/减法要求两个向量维度相同，否则返回错误
- 点积要求两个向量维度相同
- 叉积仅支持3D向量
- 混合积仅支持3D向量
- 向量与标量的加减返回错误（类型不匹配）
- 向量与向量的乘法返回错误（应使用 `dot()` 函数）

**向量算术运算**（在 `Value` 运算符重载中实现）：

| 运算 | 语法 | 实现 |
|------|------|------|
| 加法 | `VerA+VerB` | 逐分量相加，维度必须相同 |
| 减法 | `VerA-VerB` | 逐分量相减，维度必须相同 |
| 标量乘法 | `3*VerA` 或 `VerA*3` | 标量与每个分量相乘 |
| 取负 | `-VerA` | 每个分量取负 |

**向量函数**：

| 函数 | 签名 | 算法 | 约束 |
|------|------|------|------|
| `vecmod(v)` | 1个向量 | \|v\| = √(v₁²+v₂²+...+vₙ²)，调用 `eval_sqrt` | 非空向量 |
| `dot(a, b)` | 2个向量 | Σaᵢbᵢ | 同维度，非空 |
| `cross(a, b)` | 2个向量 | (a₂b₃-a₃b₂, a₃b₁-a₁b₃, a₁b₂-a₂b₁) | 仅3D |
| `scalarmul(s, v)` | 标量+向量 | s×v（逐分量乘法） | 非空向量 |
| `mixed(a, b, c)` | 3个向量 | a·(b×c) = dot(a, cross(b, c)) | 仅3D |
| `proj(a, b)` | 2个向量 | (a·b / b·b) × b | 同维度，b≠0 |
| `decompose(a, b, c)` | 3个向量（2D） | Cramer 法则解 a=αb+βc | 2D，b,c 非零且不共线 |
| `decompose(a, b, c, d)` | 4个向量（3D） | Cramer 法则解 a=αb+βc+γd | 3D，b,c,d 不共面 |

**向量分解算法**：

2D 分解（平面向量基本定理）：给定向量 a 和两个不共线的基向量 b、c，求解 a = αb + βc。

```
| a₁ |   | b₁  c₁ | | α |       det = b₁c₂ - b₂c₁
| a₂ | = | b₂  c₂ | | β |       α = (a₁c₂ - a₂c₁) / det
                                    β = (b₁a₂ - b₂a₁) / det
```

3D 分解（空间向量基本定理）：给定向量 a 和三个不共面的基向量 b、c、d，求解 a = αb + βc + γd。

使用 Cramer 法则，计算 3×3 行列式：

```
det = b₁(c₂d₃-c₃d₂) - c₁(b₂d₃-b₃d₂) + d₁(b₂c₃-b₃c₂)
```

将 a 分别替换对应列得到 det_a、det_b、det_c，则 α = det_a/det，β = det_b/det，γ = det_c/det。

**错误处理**：
- 维度不匹配："Vector dimension mismatch in addition/subtraction"
- 向量与标量运算："Cannot add/subtract vector and scalar"
- 向量×向量："Use dot(a,b) for dot product"
- 零向量基："basis vectors cannot be zero"
- 共线基向量："basis vectors are collinear"
- 共面基向量（3D）："basis vectors are coplanar"
- 叉积非3D："cross product requires 3D vectors"
- 投影到零向量："cannot project onto zero vector"

**向量变量与标量变量的分离存储**：
- 标量变量存储在 `Evaluator::variables_`（`map<string, Value>`）
- 向量变量存储在 `Evaluator::vec_variables_`（`map<string, Value>`）
- `get_variable()` 同时查找两个映射表，优先查找标量变量
- REPL 赋值时通过 `is_vec_var()` 判断使用 `set_variable()` 还是 `set_vec_variable()`

### 5.8 逗号分隔多表达式

**语法**：在一行输入中，使用逗号分隔多个表达式，依次求值并输出每个结果。

**示例**：
```
CCalc> x=1, x+1, 2/2, 1-2*5
x = 1
2
1
-9
```

**实现方式**：

REPL 在求值前调用 `split_top_level(input)` 将输入按顶层逗号拆分为多个子表达式。该函数通过括号深度计数区分顶层逗号和括号内逗号：

```
depth = 0
遍历每个字符：
  遇到 '(' → depth++
  遇到 ')' → depth--
  遇到 ',' 且 depth == 0 → 在此处拆分
```

**逗号歧义处理**：

| 场景 | 逗号层级 | 处理方式 |
|------|---------|---------|
| `x=1, x+1` | 顶层逗号 | 拆分为两个表达式 |
| `log(2, 8)` | 函数参数逗号（depth=1） | 不拆分，作为函数参数 |
| `(1,2,3)` | 向量字面量逗号（depth=1） | 不拆分，作为向量字面量 |
| `VerA=(1,2,3), vecmod(VerA)` | 混合 | 仅在顶层逗号处拆分 |

**求值顺序**：各子表达式按从左到右的顺序依次求值，前一个表达式的赋值效果对后续表达式可见（如 `x=1, x+1` 中 `x+1` 能读取到 `x=1` 的赋值）。

**错误处理**：任一子表达式出错时，输出错误信息并继续求值后续子表达式。

**矩阵字面量中的逗号处理**：`split_top_level` 增加了 `bracket_depth` 计数，方括号内的逗号不拆分：

```
depth = 0, bracket_depth = 0
遍历每个字符：
  遇到 '(' → depth++
  遇到 ')' → depth--
  遇到 '[' → bracket_depth++
  遇到 ']' → bracket_depth--
  遇到 ',' 且 depth == 0 且 bracket_depth == 0 → 在此处拆分
```

### 5.10 矩阵运算

**矩阵字面量语法**：`[1,2;3,4]`，分号分隔行，逗号分隔列。

**矩阵算术运算**（在 `Value` 运算符重载中实现）：

| 运算 | 语法 | 实现 |
|------|------|------|
| 加法 | `A + B` | 逐元素相加，维度必须相同 |
| 减法 | `A - B` | 逐元素相减，维度必须相同 |
| 矩阵乘法 | `A * B` | 标准矩阵乘法，A 的列数 = B 的行数 |
| 标量乘法 | `3 * A` 或 `A * 3` | 标量与每个元素相乘 |
| 取负 | `-A` | 每个元素取负 |

**矩阵函数**：

| 函数 | 签名 | 算法 | 约束 |
|------|------|------|------|
| `det(m)` | 方阵 | 递归余子式展开（1×1直接返回，2×2/3×3硬编码，n≥4递归） | 方阵 |
| `inv(m)` | 方阵 | 伴随矩阵法：adj(M)/det(M)，先计算所有余子式再转置 | 方阵，det≠0 |
| `eigen(m)` | 方阵 | 2×2/3×3 特征多项式 + 数值求根 | ≤4×4 |
| `trace(m)` | 方阵 | 对角线元素之和 | 方阵 |
| `transpose(m)` | 任意 | 行列互换 | — |
| `identity(n)` | 整数 | 对角线1，其余0 | 1≤n≤100 |

**行列式计算**：

- 1×1：直接返回唯一元素
- 2×2：`ad - bc`
- 3×3：Sarrus 法则（硬编码展开）
- n×4+：按第一行余子式展开，递归调用 `eval_det`

**逆矩阵计算**：

1. 计算 det(M)，若为 0 则报错 "Matrix is singular"
2. 对每个元素 (i,j) 计算余子式（删去第 i 行第 j 列的子矩阵的行列式）
3. 乘以符号 (-1)^(i+j) 得到代数余子式
4. 转置得到伴随矩阵 adj(M)
5. 每个元素除以 det(M)

**特征值计算**：

- 2×2：特征多项式 λ² - tr(A)λ + det(A) = 0，用二次方程求根
- 3×3：特征多项式 λ³ - tr(A)λ² + ((tr(A)² - tr(A²))/2)λ - det(A) = 0，用 Cardano 公式
- 使用 `double` 精度计算，结果格式化输出

### 5.11 统计函数

| 函数 | 签名 | 算法 |
|------|------|------|
| `mean(a,b,...)` | 可变参数 | Σxᵢ/n，支持向量参数展开 |
| `variance(a,b,...)` | 可变参数 | Σ(xᵢ-x̄)²/(n-1)，样本方差（Bessel 校正） |
| `stddev(a,b,...)` | 可变参数 | √variance，调用 `eval_sqrt` |
| `median(a,b,...)` | 可变参数 | 排序后取中位数，偶数个取平均 |

**实现细节**：
- `mean` 使用 Value 精确算术（BigRat），结果可能为精确有理数
- `variance` 和 `median` 使用 `double` 精度计算（通过 `val_to_double` 转换）
- `stddev` 对 variance 结果调用 `eval_sqrt`，可能返回精确根式
- 向量参数自动展开：`mean((1,2,3))` 等价于 `mean(1,2,3)`

### 5.12 自定义函数

**定义语法**：`f(x) := expr`，`g(x,y) := expr`

**解析过程**（`parse_custom_function`）：
1. 检测 `:=`（COLON_EQUAL token），将输入分为左部和右部
2. 左部解析：提取函数名和参数列表
3. 右部解析：通过 Lexer + Parser 构造 AST
4. 存入 `Evaluator::user_functions_` 映射表

**存储结构**：
```cpp
std::map<std::string, std::pair<std::vector<std::string>, ASTPtr>> user_functions_;
// 键：函数名，值：(参数名列表, 函数体AST)
```

**调用过程**（`eval_function`）：
1. 检查 `has_user_function(name)`
2. 验证参数个数匹配
3. 对函数体 AST 进行参数替换：递归遍历 AST，将每个参数变量节点替换为对应实参的数值节点
4. 对替换后的 AST 调用 `eval_node` 求值

**限制**：
- 函数名不能与内置函数名冲突
- 参数替换仅替换 VARIABLE 节点，不替换函数名中的变量
- 递归调用自定义函数不支持（会导致无限循环）

### 5.13 表达式化简

**函数**：`simplify(expr)` — 求值并显示化简结果

**实现方式**：
1. 对表达式 AST 求值得到数值结果
2. 对于除法表达式，分别求值分子和分母，显示 `分子/分母 = 结果`
3. 对于其他表达式，显示 `原式 = 结果`
4. 使用 `ast_to_string(node)` 将 AST 转为可读字符串

**示例**：
- `simplify(6/3)` → `6/3 = 2`
- `simplify(2+3)` → `2+3 = 5`

### 5.9 方程求解

**函数**：`solve(a, b, ...)` — 传入多项式各项系数（从最高次到常数项），返回所有根的向量。

| 参数个数 | 方程类型 | 示例 |
|---------|---------|------|
| 3 | 二次 ax²+bx+c=0 | `solve(1,2,1)` → `(-1, -1)` |
| 4 | 三次 ax³+bx²+cx+d=0 | `solve(1,0,0,-1)` → `(1, -0.5+0.866i, -0.5-0.866i)` |
| 5 | 四次 ax⁴+bx³+cx²+dx+e=0 | `solve(1,-10,35,-50,24)` → `(4, 1, 3, 2)` |

**二次方程（精确根式解）**：

使用 Value 精确算术系统求解，保留根式精确形式：
- 判别式 D = b² - 4ac
- √D 通过 `eval_sqrt()` 计算，返回精确的 SurdsExpr
- 根 x = (-b ± √D) / (2a)
- 当 D < 0 时，√D 返回复数，自动得到复数根

**示例**：
- `solve(1,2,1)` → `(-1, -1)` — 重根
- `solve(1,-1,-1)` → `(sqrt(5)/2+1/2, -sqrt(5)/2+1/2)` — 黄金比例，精确根式
- `solve(1,0,1)` → `(i, -i)` — 纯虚数根

**三次方程（Cardano 公式）**：

使用 depressed cubic 变换 x = t - b/(3a)，化为 t³ + pt + q = 0：
- p = (3ac - b²) / (3a²)
- q = (2b³ - 9abc + 27a²d) / (27a³)

根据判别式 Q = p³/27 + q²/4 选择算法：
- Q > 0：Cardano 公式，一个实根两个共轭复根
  - α = ∛(-q/2 + √Q), β = ∛(-q/2 - √Q)
  - x₁ = α + β, x₂ = ωα + ω̄β, x₃ = ω̄α + ωβ（ω = -1/2 + i√3/2）
- Q < 0（casus irreducibilis）：三角法
  - θ = arccos(-3q√3 / (2(-p)^(3/2))) / 3
  - xₖ = 2√(-p/3) cos(θ - 2kπ/3), k=0,1,2
- Q = 0：重根情况
  - x₁ = 2∛(-q/2), x₂ = x₃ = -∛(-q/2)

**四次方程（Ferrari 法）**：

通过 Tschirnhaus 变换 x = y - b/(4a) 化为 depressed quartic y⁴ + py² + qy + r = 0：
- p = c - 3b²/8
- q = b³/8 - bc/2 + d
- r = -3b⁴/256 + b²c/16 - bd/4 + e

当 q ≈ 0 时，化为二次方程 (y²)² + p(y²) + r = 0。

一般情况，求解三次预解方程 8s³ - 4ps² - 8rs + (4pr - q²) = 0，取实根 s，然后分解为两个二次方程：
- y² + ty + (s+u) = 0
- y² - ty + (s-u) = 0
- 其中 t = √(2s-p), u = -q/(2t)

**数值精度**：三次和四次方程使用 `double` 精度计算，结果通过 `BigFloat(double)` 构造函数转为高精度表示。二次方程使用 Value 精确算术，无精度损失。

**退化处理**：
- a=0 时二次方程退化为一次方程 bx+c=0，返回单根 -c/b
- a=b=0 时，若 c=0 返回"无穷多解"，否则返回"无解"

### 5.14 欧拉函数

**函数**：`euler_phi(n)` 或 `phi(n)` — 计算 Euler's totient function φ(n)

**定义**：φ(n) 为小于 n 且与 n 互素的正整数个数。

**算法**：基于质因数分解的公式 φ(n) = n × Π(1 - 1/p)，其中 p 遍历 n 的所有不同质因子。

**实现**：
1. `result = n`
2. 对 n 进行试除法分解
3. 对每个不同的质因子 d：`result -= result / d`
4. 若最后 `temp > 1`（剩余一个大质因子）：`result -= result / temp`

**示例**：
- `euler_phi(1)` → 1
- `euler_phi(10)` → 4（1,3,7,9 与 10 互素）
- `euler_phi(30)` → 8
- `phi(100)` → 40

### 5.15 函数值表

**函数**：`table(expr, var, from, to[, step])` — 输出函数在不同变量值下的结果

**实现**：
1. 解析参数：表达式、变量名、起始值、终止值、步长（默认 1）
2. 从 `from` 到 `to` 以 `step` 为步长遍历
3. 每个值通过 `substitute(expr, var, val)` 求值
4. 格式化输出：`var = value => result`

**示例**：
- `table(x^2, x, 0, 5)` → x=0→0, x=1→1, ..., x=5→25
- `table(sin(x), x, 0, 6, 1)` → sin(x) 在 x=0,1,...,6 处的值

### 5.16 拉格朗日乘数法

**函数**：`lagrange(f, g, x, y)` — 在约束 g(x,y)=0 下求 f(x,y) 的极值

**算法**：
1. 使用 `sdiff` 计算偏导数：∂f/∂x, ∂f/∂y, ∂g/∂x, ∂g/∂y
2. 输出方程组：∂f/∂x = λ·∂g/∂x, ∂f/∂y = λ·∂g/∂y, g(x,y) = 0
3. 在 [-3,3]×[-3,3] 网格上以步长 1 生成初始点
4. 对每个初始点进行最多 20 次 Newton 迭代：
   - 计算约束函数值 g(x,y)
   - 计算偏导数值
   - 根据 ∇f 和 ∇g 的方向调整 x, y
   - 收敛条件：|dx| < 1e-10, |dy| < 1e-10, |g| < 1e-8
5. 去重：距离 < 0.01 的点视为同一驻点
6. 输出所有驻点及其函数值

**实现细节**：
- 求值时临时降低精度到 15 位以加速计算
- 使用 `BigRat` 近似值代替 `BigFloat` 避免负数幂问题
- `eval_at` lambda 函数通过 `set_variable` + `eval_node` 求值

**示例**：
- `lagrange(x+y, x^2+y^2-1, x, y)` → 在单位圆上找 x+y 的极值
  - 最大值 f=1 在 (0,1) 和 (1,0)
  - 最小值 f=-1 在 (-1,0) 和 (0,-1)

### 5.17 符号求导化简

**化简器**：`simplify_ast(ASTPtr)` — 对 AST 进行代数化简

**化简规则**：

| 规则 | 示例 |
|------|------|
| 0 + x → x | `0+x` → `x` |
| x + 0 → x | `x+0` → `x` |
| x - 0 → x | `x-0` → `x` |
| 0 - x → -x | `0-x` → `-x` |
| 0 * x → 0 | `0*x` → `0` |
| 1 * x → x | `1*x` → `x` |
| -1 * x → -x | `-1*x` → `-x` |
| x / 1 → x | `x/1` → `x` |
| 0 / x → 0 | `0/x` → `0` |
| x ^ 0 → 1 | `x^0` → `1` |
| x ^ 1 → x | `x^1` → `x` |
| 0 ^ x → 0 | `0^x` → `0` |
| x - x → 0 | `x-x` → `0` |
| x / x → 1 | `x/x` → `1` |
| 常数折叠 | `2*3` → `6`, `2+3` → `5` |
| -(-x) → x | 双重取消除去 |
| x - (-y) → x + y | 负号吸收 |
| 系数合并 | `2*(3*x)` → `6*x` |

**AST 相等性**：`ast_equal(a, b)` 递归比较两棵 AST 是否结构完全相同，用于检测 x-x 和 x/x 模式。

**优先级感知输出**：`ast_to_string` 根据运算符优先级（+/-=1, */÷=2, ^=3, 其他=4）决定是否添加括号，避免冗余括号。

## 6. 输出格式化

### 6.1 常规格式化

`format_result()` 根据值类型选择最优显示：

| 值类型 | 显示格式 | 示例 |
|--------|---------|------|
| 整数 | 纯整数 | `42` |
| 有理数（非整数） | 分数 + 小数 | `10/3 (= 3.333...)` |
| 根式 | 精确形式 + 近似值 | `sqrt(2) ~= 1.414...` |
| 含 π | 精确形式 + 近似值 | `pi/3 ~= 1.047...` |
| 含 e | 精确形式 + 近似值 | `e ~= 2.718...` |
| 纯浮点 | 小数 | `3.14159...` |
| 复数 | a + bi | `3 + 4*i` |
| 向量 | (分量, 分量, ...) | `(1, 2, 3)` |
| 矩阵 | 美化方框输出（见 6.3） | |
| 错误 | Error: 信息 | `Error: Division by zero` |
| 字符串 | 原始文本 | `2^3 * 3^2 * 5` |

### 6.2 LaTeX 格式化

通过 `latex` 命令切换 LaTeX 输出模式，`format_latex()` 将 Value 转为 LaTeX 代码：

| 值类型 | LaTeX 格式 | 示例 |
|--------|-----------|------|
| 有理数 | `\frac{p}{q}` | `1/2` → `\frac{1}{2}` |
| 整数 | 原样 | `3` → `3` |
| 根式 | `\sqrt{n}` | `sqrt(2)` → `\sqrt{2}` |
| π | `\pi` | `pi/3` → `\frac{\pi}{3}` |
| e | `e` | `3*e` → `3e` |
| 根式系数 | `\frac{p\sqrt{n}}{q}` | `sqrt(2)/2` → `\frac{\sqrt{2}}{2}` |
| 复数 | `a + bi` | `3+i` → `3 + i` |
| 向量 | `\left(...\right)` | `(1,2)` → `\left(1, 2\right)` |
| 矩阵 | `\begin{pmatrix}...\end{pmatrix}` | `[1,2;3,4]` → `\begin{pmatrix}1 & 2 \\ 3 & 4\end{pmatrix}` |

**实现方式**：
- `Value::to_latex()` 和 `Evaluator::format_latex()` 协同工作
- 对 SurdsExpr 的每个 term 分别格式化，处理正负号和系数
- 矩阵使用 `&` 分隔列，`\\` 分隔行
- REPL 中通过 `evaluator.latex_mode()` 标志控制输出路径

### 6.3 矩阵美化输出

`format_pretty_matrix()` 使用 Unicode 方框绘制字符对齐输出矩阵：

```
┌─────┬─────┐
│  1  │  2  │
├─────┼─────┤
│  3  │  4  │
└─────┴─────┘
```

**实现方式**：
1. 计算每列最大宽度 `col_widths[]`
2. 每个单元格居中对齐（`left_pad + content + right_pad`）
3. 使用 Unicode 方框绘制字符：
   - `┌ ┐ └ ┘` 四角
   - `├ ┤` 左右 T 形
   - `─` 水平线
   - `│` 竖线
   - `┬ ┴` 上下 T 形
4. 行间用 `├─┼─┤` 分隔，首行上方 `┌─┬─┐`，末行下方 `└─┴─┘`

**触发条件**：当 `result.is_matrix()` 时自动使用美化输出，而非 `format_result`。

## 7. 性能特征

### 7.1 精度与性能

- 默认精度：110 位有效数字（`g_precision = 110`）
- 可通过 `prec N` 命令调整（1 ~ 10000）
- BigInt 运算为 O(n²) 学校算法，对 100 位精度足够高效
- BigFloat 数学函数使用硬编码常量（π, e, ln2 各 1000+ 位），首次调用后缓存

### 7.2 已知限制

- 数值微分为近似值，精度约 60~70 位有效数字
- 定积分为近似值（Simpson 法则），精度取决于区间数
- `e^2` 等超越数幂次无法给出精确符号表示
- 符号求导不做化简，输出可能包含冗余项（如 `1*0`）
- 积分表仅支持基本的模式匹配，不支持换元积分和分部积分
- 极限计算为数值逼近，对振荡函数可能不准确
- 自定义函数不支持递归调用
- 特征值计算使用 double 精度，大矩阵可能有精度损失

## 8. REPL 交互界面

### 8.1 自定义输入行（readline_custom）

**设计动机**：Windows 标准输入不支持行内编辑、历史导航和语法高亮，因此使用 `_getch()` 实现自定义输入行。

**实现方式**：
- 使用 `_getch()` 逐字符读取，不经过标准输入缓冲
- 维护 `line` 字符串和 `cursor` 位置
- 每次修改后重绘整行（`\r` + prompt + 高亮内容 + 清除残余 + 定位光标）

**按键处理**：

| 按键 | 扫描码 | 行为 |
|------|--------|------|
| Enter | `\r` / `\n` | 提交输入，加入历史 |
| Backspace | `\b` / `127` | 删除光标前字符 |
| ↑ | `0x00/0xE0` + `0x48` | 历史上翻 |
| ↓ | `0x00/0xE0` + `0x50` | 历史下翻 |
| ← | `0x00/0xE0` + `0x4B` | 光标左移 |
| → | `0x00/0xE0` + `0x4D` | 光标右移 |
| Tab | `\t` | 函数名补全 |
| 可打印字符 | `32~126` | 插入到光标位置 |

**Windows 控制台特殊处理**：
- `ENABLE_VIRTUAL_TERMINAL_PROCESSING` 标志启用 ANSI 转义序列支持
- `SetConsoleOutputCP(CP_UTF8)` / `SetConsoleCP(CP_UTF8)` 设置 UTF-8 编码

### 8.2 历史记录导航

**存储**：`std::vector<std::string> history`，`int history_idx` 指向当前位置。

**行为**：
- 每次提交非空输入时，检查是否已存在于历史中，避免重复
- `history_idx` 初始指向 `history.size()`（最新位置之后）
- ↑ 键：`history_idx--`，显示 `history[history_idx]`
- ↓ 键：`history_idx++`，若超出范围则清空输入
- 历史浏览时修改输入不会覆盖原始历史记录

### 8.3 Tab 补全

**触发**：按下 Tab 键时，提取光标前的当前单词，在 `all_function_names` 中搜索匹配项。

**算法**：
1. 从 `cursor-1` 向左扫描，收集连续的字母数字字符作为 `word`
2. 在 `all_function_names` 中查找前缀匹配 `word` 的所有函数名
3. 若唯一匹配：直接补全为 `funcname(`
4. 若多个匹配：补全到最长公共前缀
5. 无匹配：不做任何操作

**示例**：
- 输入 `si` + Tab → `sin(`
- 输入 `co` + Tab → `cos(`（公共前缀为 `cos`，因为 `comb` 也匹配但公共前缀为 `co`）
- 输入 `std` + Tab → `stddev(`

### 8.4 语法高亮

**实现**：`syntax_highlight()` 函数对输入字符串进行实时着色。

**颜色方案**（ANSI 转义码，默认主题）：

| 元素 | 颜色 | ANSI 码 | 示例 |
|------|------|---------|------|
| 数字 | 灰色 | `\033[2;37m` | `3.14`, `0xFF` |
| 函数名 | 绿色 | `\033[32m` | `sin`, `cos`, `det` |
| 运算符 | 黄色 | `\033[33m` | `+`, `-`, `*`, `/`, `^`, `%` |
| 常量 | 灰色 | `\033[2;37m` | `pi`, `e` |
| 其他 | 默认色 | — | 变量名 `x`, 括号等 |

**分词规则**：
- 数字：连续的数字和小数点，包括进制前缀（`0b`, `0o`, `0x`）
- 标识符：连续的字母、数字和下划线
- 运算符：`+`, `-`, `*`, `/`, `^`, `%`
- 其他字符：不着色

### 8.5 批量计算模式

**语法**：`ccalc -f file.txt`

**实现方式**：
1. 检测 `argc >= 3 && argv[1] == "-f"`
2. 打开文件，逐行读取
3. 空行和 `#` 开头的行跳过（注释）
4. 每行通过 `split_top_level` 拆分，逐表达式求值
5. 支持变量赋值、自定义函数定义
6. LaTeX 模式下输出 LaTeX 格式
7. 图形命令在批量模式下跳过

**限制**：
- 使用 `std::getline` 读取（非 `_getch`），因此不支持交互式功能
- 图形命令被跳过（无窗口环境）

### 8.6 多行输入

**实现**：`read_full_line()` 函数在 `readline_custom()` 之上实现续行逻辑。

**续行条件**：
1. **显式续行**：行尾为 `\` 时，移除 `\`，继续读取下一行（提示符变为 `...> `）
2. **自动续行**：括号不匹配时（`count_unmatched() > 0`），自动继续读取下一行

**括号计数**：`count_unmatched(s)` 统计未匹配的 `(`、`)`、`[`、`]`，返回净未匹配数。

**历史记录处理**：
- 续行中间的行不加入历史
- 最终完整输入作为一条历史记录

**示例**：
```
CCalc> solve(1, \
...> 2, 1)
(-1, -1)

CCalc> det([1,2;
...> 3,4])
-2
```

### 8.7 配色主题

**数据结构**：`ColorTheme` 结构体包含主题名称和四种颜色码（数字、函数、运算符、常量）。

**内置主题**：

| 主题名 | 数字 | 函数 | 运算符 | 常量 |
|--------|------|------|--------|------|
| `default` | dim white (灰) | green (绿) | yellow (黄) | dim white (灰) |
| `dark` | white (白) | cyan (青) | magenta (紫) | white (白) |
| `light` | black (黑) | blue (蓝) | red (红) | black (黑) |
| `neon` | bright cyan (亮青) | bright green (亮绿) | bright yellow (亮黄) | bright cyan (亮青) |
| `none` | 无色 | 无色 | 无色 | 无色 |

**切换命令**：
- `theme` — 列出可用主题和当前主题
- `theme dark` — 切换到 dark 主题

**实现**：`syntax_highlight()` 函数读取 `themes[g_theme_idx]` 获取当前颜色码。

### 8.8 中文模式

**切换命令**：`zh` 切换到中文，`en` 切换到英文。

**翻译机制**：`T(const char* en, const char* zh)` 函数根据 `g_chinese` 标志返回对应语言的字符串。

**翻译范围**：
- 欢迎信息
- 帮助文本（完整中文版 `print_help()`）
- 错误消息
- 状态消息（角度模式、精度、进制、LaTeX 模式等）
- 退出消息

**限制**：计算结果本身不翻译（数学符号通用）。

### 8.9 配置文件

**路径**：
- Windows：`%APPDATA%\ccalcconfig`
- Linux/macOS：`~/.ccalcconfig`

**格式**：`key=value`，`#` 开头为注释。

**保存的配置项**：

| 键 | 值 | 示例 |
|----|-----|------|
| `precision` | 整数 | `precision=50` |
| `angle_mode` | `deg` 或 `rad` | `angle_mode=rad` |
| `latex_mode` | `on` 或 `off` | `latex_mode=off` |
| `language` | `zh` 或 `en` | `language=zh` |
| `theme` | 主题名 | `theme=neon` |
| `output_base` | 整数 | `output_base=10` |
| `custom_func` | 函数定义 | `custom_func=f(x):=x^2+1` |

**加载时机**：程序启动时调用 `load_config()`，若文件不存在则使用默认值。

**保存时机**：输入 `quit` 或 `exit` 时调用 `save_config()`。

**自定义函数序列化**：通过 `Evaluator::get_user_func_string(name)` 获取函数定义字符串，保存为 `custom_func=名(参数):=表达式` 格式。加载时通过 `parse_custom_function()` 重新解析和注册。

## 9. 编译与构建

**依赖**：仅 C++17 标准库，无外部依赖

**编译器要求**：支持 C++17 的 GCC/Clang/MSVC

**构建命令**：
```bash
make        # 编译
make clean  # 清理
```

**编译选项**：`-std=c++17 -O2 -Wall -Wextra`

**Windows 特殊处理**：
- `SetConsoleOutputCP(CP_UTF8)` / `SetConsoleCP(CP_UTF8)` 设置控制台 UTF-8 编码
- `std::flush` 确保提示符即时显示

## 10. 函数图像绘制

### 10.1 概述

CCalc 提供独立的函数图像绘制窗口 `ccalc_graph.exe`，基于 Dear ImGui + DirectX 11 实现。支持四种绘图类型：

| 类型 | REPL 语法 | 命令行 | 示例 |
|------|-----------|--------|------|
| 显式 y=f(x) | `graph(sin(x))` | `ccalc_graph "sin(x)"` | 正弦曲线 |
| 隐函数 f(x,y)=0 | `graph_implicit(x^2+y^2-1)` | `ccalc_graph -i "x^2+y^2-1"` | 单位圆 |
| 参数方程 | `graph_param(t, cos(t), sin(t))` | `ccalc_graph -p "cos(t)" "sin(t)"` | 圆 |
| 极坐标 r=f(θ) | `graph_polar(1+cos(theta))` | `ccalc_graph -l "1+cos(theta)"` | 心形线 |

### 10.2 架构设计

```
graph.cpp（独立可执行文件）
├── 命令行参数解析（argc/argv → init_funcs）
│   ├── 默认：显式 y=f(x)
│   ├── -i：隐函数 f(x,y)=0
│   ├── -p x_expr y_expr：参数方程
│   └── -l：极坐标 r=f(theta)
├── Win32 窗口创建 + DirectX 11 渲染后端
├── ImGui 界面
│   ├── 设置面板（默认隐藏，通过 "Settings >>" 按钮切换）
│   │   ├── 函数列表（类型选择 + 输入框 + 可见性开关 + 颜色标识）
│   │   ├── 参数范围（param/polar 的 t_min, t_max）
│   │   ├── 视图范围控制
│   │   ├── 显示选项（网格、坐标）
│   │   └── 角度模式
│   ├── 右侧绘图画布
│   │   ├── 坐标轴和网格
│   │   ├── 四种绘图算法分发
│   │   ├── 鼠标悬停坐标显示
│   │   └── 交互控制（缩放、平移）
│   └── 浮动 Settings 按钮（overlay 窗口）
└── fast_eval() — 基于 double 的快速求值器
    └── EvalVars { x, y, t, theta } 支持四种变量
```

### 10.3 快速求值器（fast_eval）

**设计动机**：CCalc 的核心求值器使用高精度 BigFloat 运算，每次求值涉及字符串构造、高精度计算和格式转换，对实时绘图（每帧需 ~1000 次求值）来说太慢。

**解决方案**：实现一个基于原生 `double` 的快速 AST 遍历求值器，完全绕过高精度运算系统。

**实现方式**：
- 递归遍历 AST 节点
- 每个节点类型直接使用 `double` 运算
- 支持的节点类型：
  - `NUMBER`：通过 `BigRat::to_string()` 转为 double（仅在 AST 缓存失效时执行一次）
  - `CONSTANT`：`pi` → `M_PI`, `e` → `M_E`
  - `VARIABLE`：`x`, `y`, `t`, `theta` → EvalVars 对应字段
  - `BINOP`：`+`, `-`, `*`, `/`, `^`, `%`
  - `UNARYOP`：一元负号
  - `FACTORIAL`：迭代计算（限制 ≤170）
  - `FUNCTION`：30+ 个数学函数的 double 实现

**角度模式处理**：`fast_eval` 接受 `deg_mode` 参数，三角函数在度数模式下自动进行弧度转换。

**不支持的函数**：`factor`, `gcd`, `lcm`, `P`, `C`, `convert`, `solve` 等返回 NAN（不适用于绘图场景）。

### 10.4 AST 缓存机制

**问题**：如果每次采样都重新解析表达式字符串（Lexer + Parser），每帧需解析数千次，导致窗口卡死。

**解决方案**：`FuncEntry` 结构缓存解析结果：
- `cached_ast: ASTPtr`：缓存的 AST 指针
- `cached_expr_str: string`：缓存对应的表达式字符串
- `ensure_ast()`：仅当表达式变化时重新解析

**缓存失效**：当用户修改输入框中的表达式时，`ensure_ast()` 检测到 `expr != cached_expr_str`，触发重新解析。

### 10.5 绘图算法

#### 10.5.1 显式函数 y=f(x)

**采样策略**：
- 采样点数 = 画布像素宽度（通常 700~1500 点）
- 范围：[x_min, x_max]，均匀分布
- 每个采样点调用 `fast_eval(ast, EvalVars{x}, deg_mode)`

**不连续性检测**：
- 当 `|y - prev_y| > (y_max - y_min) * 10` 时，认为函数不连续（如 `tan(x)` 的渐近线）
- 不连续处断开折线，避免竖直连线

#### 10.5.2 隐函数 f(x,y)=0（Marching Squares）

**算法**：Marching Squares 等值线提取
1. 在视图范围内建立 N×M 网格（约 100~300 点/轴）
2. 在每个网格点计算 f(x,y) 值
3. 对每个网格单元，根据四角值的正负确定 4-bit 案例码
4. 案例码 0 和 15（全正/全负）跳过
5. 案例码有 2 条边穿越时，线性插值求穿越点，画一条线段
6. 案例码有 4 条边穿越时（鞍点），用四角平均值判断连接方式
7. 线性插值公式：`t = f0 / (f0 - f1)`，位置 = p0 + t * (p1 - p0)

**优势**：能绘制任意隐函数曲线，如圆、椭圆、双曲线等

#### 10.5.3 参数方程 (x(t), y(t))

**采样策略**：
- 参数 t 在 [t_min, t_max] 范围内均匀采样 2000 个点
- 默认范围：[0, 2π]
- 每个采样点计算 x(t) 和 y(t)，直接映射到屏幕坐标
- 跳跃检测：相邻屏幕点距离超过 800 像素时断开折线

#### 10.5.4 极坐标 r=f(θ)

**采样策略**：
- 参数 θ 在 [t_min, t_max] 范围内均匀采样 2000 个点
- 默认范围：[0, 2π]
- 每个采样点计算 r(θ)，转换为笛卡尔坐标 x=r·cos(θ), y=r·sin(θ)
- 跳跃检测同参数方程

**自动 Y 范围**：
- 显式函数：200 个采样点预扫描
- 参数方程/极坐标：500 个采样点预扫描
- 取 Y 值范围的 10% 边距（最小 0.5）
- 在绘制网格和坐标轴之前计算

**网格间距**：使用"nice number"算法选择网格间距：
```
rough = range / 10
选择 1, 2, 5, 10 中最接近 rough 的数量级
```

### 10.6 交互功能

| 操作 | 实现 |
|------|------|
| 滚轮缩放 | 以鼠标位置为中心，仅缩放 X 轴，Y 轴由 auto_y 调整 |
| 左键拖拽平移 | 记录拖拽起始位置和范围，按像素偏移量换算 |
| 鼠标悬停 | 显示十字准线 + 坐标值 + 各显式函数在该 x 处的值和标记点 |
| Settings 按钮 | 左上角浮动 overlay 按钮，控制设置面板的显隐 |

### 10.7 REPL 集成

**语法**：
- `graph(expr)` — 显式函数 y=expr
- `graph_implicit(expr)` — 隐函数 f(x,y)=0
- `graph_param(t, x_expr, y_expr)` — 参数方程
- `graph_polar(r_expr)` — 极坐标

**实现方式**：
- `parse_graph_cmd()` 函数在 REPL 层面拦截 `graph(...)` / `graph_implicit(...)` / `graph_param(...)` / `graph_polar(...)` 输入
- 提取括号内的原始文本（不经过求值），构造命令行参数
- `graph_param` 使用 `split_args()` 在括号内按顶层逗号分割参数
- 使用 `start` 命令异步启动绘图窗口，REPL 不阻塞

**命令行参数**：
- `ccalc_graph.exe "sin(x)"` — 直接绘制 sin(x)
- `ccalc_graph.exe` — 无参数时默认绘制 sin(x)
- 窗口标题显示 `y = expr`

### 10.8 构建系统

**Makefile 新增目标**：
- `ccalc_graph.exe`：链接 `graph.o` + 计算器核心对象 + ImGui 对象 + DirectX 库
- ImGui 编译：自动编译 `imgui-1.92.7/` 下的 7 个源文件（imgui.cpp, imgui_draw.cpp, imgui_tables.cpp, imgui_widgets.cpp, imgui_demo.cpp, imgui_impl_win32.cpp, imgui_impl_dx11.cpp）
- 链接库：`-ld3d11 -ld3dcompiler -ldxgi -ldwmapi -lgdi32 -luser32 -lshell32`

**依赖关系**：
- `graph.cpp` 依赖 `ccalc.h`（AST 定义、Lexer、Parser）
- 不依赖 `evaluator.cpp`（使用自己的 `fast_eval`）
- 依赖 ImGui 头文件和 Win32/DX11 头文件
