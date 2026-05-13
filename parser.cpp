#include "ccalc.h"

static const std::vector<std::string> known_functions = {
    "sin", "cos", "tan", "asin", "acos", "atan",
    "sinh", "cosh", "tanh",
    "sqrt", "cbrt", "nrt",
    "abs", "ln", "lg", "log",
    "fact", "perm", "P", "comb", "C",
    "gcd", "lcm", "factor", "euler_phi", "euler", "phi",
    "int", "diff", "sum", "prod",
    "complex", "re", "im", "conj", "arg", "mod",
    "rand", "randint",
    "convert", "solve",
    "deg", "rad",
    "vecmod", "dot", "cross", "scalarmul", "mixed", "proj", "decompose",
    "matrix", "det", "inv", "eigen", "trace", "transpose", "identity",
    "mean", "stddev", "variance", "median",
    "simplify",
    "taylor", "limit", "inttable", "recur", "table", "lagrange"
};

static bool is_known_function(const std::string& name) {
    for (auto& f : known_functions) {
        if (f == name) return true;
    }
    return false;
}

char Lexer::peek() const {
    if (pos_ >= (int)input_.size()) return '\0';
    return input_[pos_];
}

char Lexer::advance() {
    if (pos_ >= (int)input_.size()) return '\0';
    return input_[pos_++];
}

void Lexer::skip_whitespace() {
    while (pos_ < (int)input_.size() && std::isspace(input_[pos_]))
        pos_++;
}

Token Lexer::read_number() {
    std::string num;
    int start = pos_;
    bool has_dot = false;
    if (pos_ < (int)input_.size() && input_[pos_] == '0' &&
        pos_ + 1 < (int)input_.size()) {
        char next = input_[pos_ + 1];
        if (next == 'b' || next == 'B') {
            pos_ += 2;
            while (pos_ < (int)input_.size() && 
                   (input_[pos_] == '0' || input_[pos_] == '1')) {
                num += input_[pos_]; pos_++;
            }
            BigInt val = BigInt::from_base_string(num, 2);
            return Token(TokenType::NUMBER, val.to_string(), start);
        }
        if (next == 'o' || next == 'O') {
            pos_ += 2;
            while (pos_ < (int)input_.size() && input_[pos_] >= '0' && input_[pos_] <= '7') {
                num += input_[pos_]; pos_++;
            }
            BigInt val = BigInt::from_base_string(num, 8);
            return Token(TokenType::NUMBER, val.to_string(), start);
        }
        if (next == 'x' || next == 'X') {
            pos_ += 2;
            while (pos_ < (int)input_.size() && 
                   std::isxdigit(input_[pos_])) {
                num += input_[pos_]; pos_++;
            }
            BigInt val = BigInt::from_base_string(num, 16);
            return Token(TokenType::NUMBER, val.to_string(), start);
        }
    }
    while (pos_ < (int)input_.size()) {
        char c = input_[pos_];
        if (c == '.' && !has_dot) {
            if (pos_ + 1 < (int)input_.size() && std::isdigit(input_[pos_ + 1])) {
                has_dot = true;
                num += c;
                pos_++;
            } else {
                break;
            }
        } else if (std::isdigit(c)) {
            num += c;
            pos_++;
        } else {
            break;
        }
    }
    return Token(TokenType::NUMBER, num, start);
}

Token Lexer::read_identifier() {
    std::string id;
    int start = pos_;
    while (pos_ < (int)input_.size() &&
           (std::isalnum(input_[pos_]) || input_[pos_] == '_')) {
        id += input_[pos_];
        pos_++;
    }
    return Token(TokenType::IDENTIFIER, id, start);
}

std::vector<Token> Lexer::tokenize() {
    std::vector<Token> tokens;
    while (pos_ < (int)input_.size()) {
        skip_whitespace();
        if (pos_ >= (int)input_.size()) break;
        char c = input_[pos_];
        if (std::isdigit(c)) {
            tokens.push_back(read_number());
        } else if (std::isalpha(c) || c == '_') {
            tokens.push_back(read_identifier());
        } else {
            switch (c) {
            case '+': tokens.push_back(Token(TokenType::PLUS, "+", pos_)); pos_++; break;
            case '-': tokens.push_back(Token(TokenType::MINUS, "-", pos_)); pos_++; break;
            case '*': tokens.push_back(Token(TokenType::STAR, "*", pos_)); pos_++; break;
            case '/': tokens.push_back(Token(TokenType::SLASH, "/", pos_)); pos_++; break;
            case '^': tokens.push_back(Token(TokenType::CARET, "^", pos_)); pos_++; break;
            case '%': tokens.push_back(Token(TokenType::PERCENT, "%", pos_)); pos_++; break;
            case '(': tokens.push_back(Token(TokenType::LPAREN, "(", pos_)); pos_++; break;
            case ')': tokens.push_back(Token(TokenType::RPAREN, ")", pos_)); pos_++; break;
            case '[': tokens.push_back(Token(TokenType::LBRACKET, "[", pos_)); pos_++; break;
            case ']': tokens.push_back(Token(TokenType::RBRACKET, "]", pos_)); pos_++; break;
            case ',': tokens.push_back(Token(TokenType::COMMA, ",", pos_)); pos_++; break;
            case ';': tokens.push_back(Token(TokenType::SEMICOLON, ";", pos_)); pos_++; break;
            case '!': tokens.push_back(Token(TokenType::BANG, "!", pos_)); pos_++; break;
            case ':':
                if (pos_ + 1 < (int)input_.size() && input_[pos_ + 1] == '=') {
                    tokens.push_back(Token(TokenType::COLON_EQUAL, ":=", pos_)); pos_ += 2;
                } else {
                    tokens.push_back(Token(TokenType::ERROR, ":", pos_)); pos_++;
                }
                break;
            case '=': tokens.push_back(Token(TokenType::EQUAL, "=", pos_)); pos_++; break;
            case '<':
                if (pos_ + 1 < (int)input_.size() && input_[pos_ + 1] == '=') {
                    tokens.push_back(Token(TokenType::LE, "<=", pos_)); pos_ += 2;
                } else {
                    tokens.push_back(Token(TokenType::LT, "<", pos_)); pos_++;
                }
                break;
            case '>':
                if (pos_ + 1 < (int)input_.size() && input_[pos_ + 1] == '=') {
                    tokens.push_back(Token(TokenType::GE, ">=", pos_)); pos_ += 2;
                } else {
                    tokens.push_back(Token(TokenType::GT, ">", pos_)); pos_++;
                }
                break;
            default:
                tokens.push_back(Token(TokenType::ERROR, std::string(1, c), pos_));
                pos_++;
                break;
            }
        }
    }
    tokens.push_back(Token(TokenType::END_OF_INPUT, "", pos_));
    std::vector<Token> result;
    for (size_t i = 0; i < tokens.size(); i++) {
        result.push_back(tokens[i]);
        if (i + 1 < tokens.size()) {
            auto& cur = tokens[i];
            auto& next = tokens[i + 1];
            bool need_mul = false;
            if (cur.type == TokenType::NUMBER && next.type == TokenType::IDENTIFIER)
                need_mul = true;
            if (cur.type == TokenType::NUMBER && next.type == TokenType::LPAREN)
                need_mul = true;
            if (cur.type == TokenType::RPAREN && next.type == TokenType::LPAREN)
                need_mul = true;
            if (cur.type == TokenType::RPAREN && next.type == TokenType::NUMBER)
                need_mul = true;
            if (cur.type == TokenType::RPAREN && next.type == TokenType::IDENTIFIER)
                need_mul = true;
            if (cur.type == TokenType::NUMBER && next.type == TokenType::IDENTIFIER) {
                if (next.text == "i" || next.text == "I") {
                    need_mul = false;
                }
            }
            if (need_mul) {
                result.push_back(Token(TokenType::STAR, "*", -1));
            }
        }
    }
    return result;
}

const Token& Parser::peek() const {
    if (pos_ >= tokens_.size()) {
        static Token end(TokenType::END_OF_INPUT, "", 0);
        return end;
    }
    return tokens_[pos_];
}

Token Parser::advance() {
    if (pos_ >= tokens_.size()) return Token(TokenType::END_OF_INPUT, "", 0);
    return tokens_[pos_++];
}

bool Parser::match(TokenType t) {
    if (check(t)) {
        advance();
        return true;
    }
    return false;
}

bool Parser::check(TokenType t) const {
    return peek().type == t;
}

ASTPtr Parser::parse() {
    auto result = expression();
    if (!check(TokenType::END_OF_INPUT)) {
        throw CalcError("Unexpected token: " + peek().text);
    }
    return result;
}

ASTPtr Parser::expression() {
    return comparison();
}

ASTPtr Parser::comparison() {
    auto left = additive();
    while (check(TokenType::LT) || check(TokenType::GT) ||
           check(TokenType::LE) || check(TokenType::GE) ||
           check(TokenType::EQUAL) || check(TokenType::NEQ)) {
        auto tok = advance();
        auto right = additive();
        char op = '=';
        if (tok.type == TokenType::LT) op = '<';
        else if (tok.type == TokenType::GT) op = '>';
        else if (tok.type == TokenType::LE) op = 'L';
        else if (tok.type == TokenType::GE) op = 'G';
        else if (tok.type == TokenType::EQUAL) op = '=';
        left = ASTNode::make_binop(op, left, right);
    }
    return left;
}

ASTPtr Parser::additive() {
    auto left = multiplicative();
    while (check(TokenType::PLUS) || check(TokenType::MINUS)) {
        auto tok = advance();
        auto right = multiplicative();
        left = ASTNode::make_binop(tok.text[0], left, right);
    }
    return left;
}

ASTPtr Parser::multiplicative() {
    auto left = power();
    while (check(TokenType::STAR) || check(TokenType::SLASH) || check(TokenType::PERCENT)) {
        auto tok = advance();
        auto right = power();
        left = ASTNode::make_binop(tok.text[0], left, right);
    }
    return left;
}

ASTPtr Parser::power() {
    auto left = unary();
    if (check(TokenType::CARET)) {
        advance();
        auto right = power();
        return ASTNode::make_binop('^', left, right);
    }
    return left;
}

ASTPtr Parser::unary() {
    if (check(TokenType::MINUS)) {
        advance();
        auto operand = unary();
        return ASTNode::make_unaryop('-', operand);
    }
    if (check(TokenType::PLUS)) {
        advance();
        return unary();
    }
    return postfix();
}

ASTPtr Parser::postfix() {
    auto left = primary();
    while (check(TokenType::BANG)) {
        advance();
        left = ASTNode::make_factorial(left);
    }
    return left;
}

ASTPtr Parser::primary() {
    if (check(TokenType::NUMBER)) {
        auto tok = advance();
        return ASTNode::make_num(BigRat(tok.text));
    }
    if (check(TokenType::IDENTIFIER)) {
        auto tok = advance();
        std::string name = tok.text;
        if (name == "pi" || name == "PI" || name == "Pi") {
            return ASTNode::make_const("pi");
        }
        if (name == "e" || name == "E") {
            if (check(TokenType::IDENTIFIER)) {
                auto next = peek();
                if (next.text.size() > 1 || !std::isdigit(next.text[0])) {
                    return ASTNode::make_const("e");
                }
            } else {
                return ASTNode::make_const("e");
            }
        }
        if (check(TokenType::LPAREN) || is_known_function(name)) {
            if (check(TokenType::LPAREN)) {
                return function_call(name);
            }
            if (name == "ans") return ASTNode::make_var("ans");
            return ASTNode::make_var(name);
        }
        return ASTNode::make_var(name);
    }
    if (check(TokenType::LPAREN)) {
        advance();
        auto first = expression();
        if (check(TokenType::COMMA)) {
            std::vector<ASTPtr> components;
            components.push_back(first);
            while (check(TokenType::COMMA)) {
                advance();
                components.push_back(expression());
            }
            if (!match(TokenType::RPAREN)) {
                throw CalcError("Expected ')' in vector literal");
            }
            return ASTNode::make_vec_literal(components);
        }
        if (!match(TokenType::RPAREN)) {
            throw CalcError("Expected ')'");
        }
        return first;
    }
    if (check(TokenType::LBRACKET)) {
        advance();
        std::vector<std::vector<ASTPtr>> rows;
        if (check(TokenType::RBRACKET)) {
            advance();
            return ASTNode::make_mat_literal(rows);
        }
        do {
            std::vector<ASTPtr> row;
            row.push_back(expression());
            while (check(TokenType::COMMA)) {
                advance();
                row.push_back(expression());
            }
            rows.push_back(row);
            if (check(TokenType::SEMICOLON)) {
                advance();
            } else {
                break;
            }
        } while (true);
        if (!match(TokenType::RBRACKET)) {
            throw CalcError("Expected ']'");
        }
        return ASTNode::make_mat_literal(rows);
    }
    throw CalcError("Unexpected token: " + peek().text);
}

ASTPtr Parser::function_call(const std::string& name) {
    advance();
    std::vector<ASTPtr> args;
    if (!check(TokenType::RPAREN)) {
        args.push_back(expression());
        while (check(TokenType::COMMA)) {
            advance();
            args.push_back(expression());
        }
    }
    if (!match(TokenType::RPAREN)) {
        throw CalcError("Expected ')' in function call");
    }
    return ASTNode::make_func(name, args);
}
