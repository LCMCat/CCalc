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
└── Makefile         # 构建脚本
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
enum Type { SURDS, FLOAT, COMPLEX, VECTOR, STRING, ERROR };
```

| 类型 | 用途 | 存储字段 |
|------|------|---------|
| SURDS | 精确值（有理数/根式/π/e） | `surds: SurdsExpr` |
| FLOAT | 浮点近似值 | `float_val: BigFloat` |
| COMPLEX | 复数 | `complex: ComplexVal` |
| VECTOR | 向量 | `vec: vector<Value>` |
| STRING | 字符串（因数分解结果） | `error_msg: string` |
| ERROR | 错误信息 | `error_msg: string` |

**类型提升规则**：
- SURDS ⊕ FLOAT → FLOAT
- SURDS/FLOAT ⊕ COMPLEX → COMPLEX
- VECTOR ⊕ VECTOR → VECTOR（同维度加减、标量乘法）
- 任何 ⊕ ERROR → ERROR

**ComplexVal**：`shared_ptr<Value> real, imag`，支持任意精度的实部和虚部。

## 4. 解析器设计

### 4.1 Lexer — 词法分析器

**Token 类型**：

```
NUMBER, IDENTIFIER, PLUS, MINUS, STAR, SLASH, CARET, PERCENT,
LPAREN, RPAREN, COMMA, BANG, EQUAL, LT, GT, LE, GE, NEQ,
END_OF_INPUT, ERROR
```

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
NUMBER     — 数值字面量 (BigRat)
CONSTANT   — 常量 (pi, e)
VARIABLE   — 变量 (x, k, ans, VerA, VerB, VerC, VerD)
BINOP      — 二元运算 (+, -, *, /, ^, %)
UNARYOP    — 一元运算 (-)
FUNCTION   — 函数调用 (sin, cos, ...)
FACTORIAL  — 阶乘后缀 (n!)
VEC_LITERAL — 向量字面量 ((1,2,3))
```

**已知函数列表**：50+ 个内置函数名，用于区分函数调用和变量名乘法。

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

**求和/连乘** — 逐项求值：
- `sum(expr, var, from, to)` — 从 from 到 to 逐整数求值并累加
- `prod(expr, var, from, to)` — 从 from 到 to 逐整数求值并累乘

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

### 5.6 进制输出

**输出进制**：通过 `base N` 命令设置（2~36），默认为 10。

**实现方式**：
- `BigInt::to_base_string(base)` 将整数转为指定进制字符串
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

## 6. 输出格式化

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
| 错误 | Error: 信息 | `Error: Division by zero` |
| 字符串 | 原始文本 | `2^3 * 3^2 * 5` |

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

## 8. 编译与构建

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
