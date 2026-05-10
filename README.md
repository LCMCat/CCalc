# CCalc

A high-precision command-line advanced calculator written in C++17, with zero external dependencies.

All core arithmetic (BigInt, BigRat, BigFloat) is implemented from scratch, supporting **100+ digit** signed numbers and **exact symbolic computation** (e.g. `sin(pi/2) = 1`, `ln(e) = 1`).

## Features

### Arithmetic
- `+  -  *  /  ^  %` with arbitrary precision
- Exact rational arithmetic: `10/3` → `10/3 (= 3.333...)`
- Implicit multiplication: `2pi`, `3(1+2)`, `(2)pi`
- Comma-separated expressions: `x=1, x+1, 2/2` → evaluates and outputs each result

### Constants
- `pi` — π, with 1000+ digit precision
- `e` — Euler's number, with 1000+ digit precision

### Variables
- `x`, `y`, `A`, `B`, `C`, `D` — user-assignable variables (default: 0)
- Assignment: `x=1`, `x=sin(pi/2)`, `A=10`
- Usage in expressions: `sin(x/2)` when `x=pi` → `1`

### Vectors
- `VerA`, `VerB`, `VerC`, `VerD` — user-assignable vector variables (default: empty)
- Assignment: `VerA=(1,2,3)`, `VerB=(1,2)`
- Vector literals: `(1,2,3)` — 3D vector, `(3,4)` — 2D vector
- Arithmetic: `VerA+VerB`, `VerA-VerB`, `3*VerA`, `-VerA` (same dimension required for ±)
- Operations:
  - `vecmod(v)` — magnitude: `vecmod((3,4))` → `5`
  - `dot(a, b)` — dot product (same dimension required)
  - `cross(a, b)` — cross product (3D only)
  - `scalarmul(s, v)` — scalar multiplication (also `s*v` or `v*s`)
  - `mixed(a, b, c)` — scalar triple product a·(b×c) (3D only)
  - `proj(a, b)` — projection of a onto b
  - `decompose(a, b, c)` — 2D decomposition: a = αb + βc, returns (α, β)
  - `decompose(a, b, c, d)` — 3D decomposition: a = αb + βc + γd, returns (α, β, γ)

### Number Bases
- Input prefixes: `0b1010` (binary), `0o77` (octal), `0xFF` (hex)
- Output base: `base 2` / `base 8` / `base 16` / `base N` (2–36)
- Example: `base 16; 255` → `FF`

### Trigonometric
- `sin`, `cos`, `tan` — with exact values at special angles
- `asin`, `acos`, `atan` — inverse trigonometric
- `sinh`, `cosh`, `tanh` — hyperbolic
- Angle modes: `deg` / `rad`

### Roots & Powers
- `sqrt(x)`, `cbrt(x)`, `nrt(n, x)` — square/cube/nth root
- `x^n` — power (exact for integer exponents)
- `n!` or `fact(n)` — factorial

### Logarithmic
- `ln(x)` — natural logarithm (`ln(e) = 1`)
- `lg(x)` — common logarithm (`lg(100) = 2`)
- `log(base, x)` — arbitrary base (`log(2, 8) = 3`)

### Number Theory
- `factor(n)` — prime factorization (`factor(360)` → `2^3 * 3^2 * 5`)
- `gcd(a, b)`, `lcm(a, b)`
- `P(n, k)` / `perm(n, k)` — permutation
- `C(n, k)` / `comb(n, k)` — combination

### Calculus
- `int(expr, var, a, b)` — definite integral (Simpson's rule)
- `diff(expr, var, point)` — numerical derivative (5-point central difference)
- `sum(expr, var, from, to)` — summation
- `prod(expr, var, from, to)` — product

### Complex Numbers
- `complex(a, b)` — construct `a + bi`
- `re(z)`, `im(z)`, `conj(z)`, `arg(z)`, `mod(z)`
- Automatic complexification: `sqrt(-1)` → `i`, `ln(-1)` → `pi*i`

### Other
- `abs(x)`, `floor(x)`, `ceil(x)`, `round(x)`, `sign(x)`
- `max(a, b)`, `min(a, b)`
- `rand()` — random in [0, 1)
- `randint(a, b)` — random integer
- `convert(val, from, to)` — unit conversion (length, mass, time, area, volume, speed, data, temperature)

## Exact Value Examples

| Input | Output |
|-------|--------|
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

## Build

### Prerequisites
- C++17 compatible compiler (GCC ≥ 9, Clang ≥ 10, MSVC ≥ 2019)
- Make

### Compile

```bash
make
```

### Clean

```bash
make clean
```

## Usage

Run the calculator:

```bash
./ccalc
```

### Interactive Session

```
CCalc v1.0 - Command-line Advanced Calculator
Type 'help' for help, 'quit' to exit.

CCalc> sin(pi/2)
1
CCalc> 2^100
1267650600228229401496703205376
CCalc> 100!
93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000
CCalc> log(2, sqrt(sin(pi/(2+5)^0)))
Error: log(0) is undefined
CCalc> complex(3,4) * complex(1,2)
-5 + 10*i
CCalc> x=pi
x = pi ~= 3.14159...
CCalc> sin(x/2)
1
CCalc> 0xFF
255
CCalc> base 16
Output base set to 16
CCalc> 255
FF
CCalc> VerA=(1,2,3)
VerA = (1, 2, 3)
CCalc> vecmod(VerA)
sqrt(14) ~= 3.74165738677394138558374873231654930176
CCalc> dot(VerA,(4,5,6))
32
CCalc> cross(VerA,(4,5,6))
(-3, 6, -3)
CCalc> 3*VerA
(3, 6, 9)
CCalc> proj((1,2),(3,4))
(33/25 (= 1.32), 44/25 (= 1.76))
CCalc> x=1, x+1, 2/2, 1-2*5
x = 1
2
1
-9
CCalc> sin(pi/2), cos(pi/3), log(2,8)
1
1/2 (= 0.5)
3
CCalc> quit
Goodbye!
```

### Commands

| Command | Description |
|---------|-------------|
| `deg` | Switch to degree mode |
| `rad` | Switch to radian mode |
| `base N` | Set output base (2–36) |
| `prec N` | Set precision to N digits (1–10000) |
| `ans` | Last answer |
| `help` | Show help |
| `quit` / `exit` | Exit |

## Architecture

```
┌──────────┐    ┌─────────┐    ┌────────┐    ┌──────────┐    ┌───────┐
│  Input   │───>│  Lexer  │───>│ Parser │───>│ Evaluator│───>│ Value │───> Output
│ (string) │    │ (tokens)│    │  (AST) │    │  (eval)  │    │(result)│
└──────────┘    └─────────┘    └────────┘    └──────────┘    └───────┘
```

### Core Types

| Type | Description |
|------|-------------|
| `BigInt` | Arbitrary-precision integer (base 10^9) |
| `BigRat` | Arbitrary-precision rational (auto-simplified) |
| `BigFloat` | Arbitrary-precision float (scientific notation) |
| `SurdsExpr` | Exact symbolic expression (radicals + π + e) |
| `Value` | Discriminated union: SURDS / FLOAT / COMPLEX / VECTOR / STRING / ERROR |

### Key Design Decisions

- **SurdsExpr** enables exact symbolic computation. Special radicand values: `-1` = π, `-2` = e, `≥2` = √n.
- **Exact-first evaluation**: operations stay in the exact (SurdsExpr) domain as long as possible, falling back to BigFloat only when necessary.
- **No external dependencies**: all high-precision arithmetic is implemented from scratch.
- **Hard-coded constants**: π, e, ln(2) are stored with 1000+ digits for maximum accuracy.

## License

MIT
