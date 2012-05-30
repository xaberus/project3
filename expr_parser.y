%token_type { token_t * }
//%token_destructor { ast_delete_token(ast, $$); }

%name "expr_parser"

%include {
  #include <assert.h>
  #include "expr.h"
  #define DEBUG
}

%token_prefix TOK_

%extra_argument { ast_t * ast }

%syntax_error {
  fprintf(stderr, "syntax error: unexpected token: '%s'\n", yyTokenName[yymajor]);
}

%fallback  OPEN LPAREN .

block ::= definition_list .

definition_list ::= .
definition_list ::= definition_list variable_definition(var) . {
  ast_deflist_append_var(ast, var);
}
definition_list ::= definition_list function_definition(func) . {
  ast_deflist_append_func(ast, func);
}

%type variable_definition { vardef_t * }
variable_definition(R) ::= NAME(name) ASSIGN expression(expr) SEMICOLON . {
  R = ast_define_variable(ast, name, expr);
}

%type function_definition { funcdef_t * }
function_definition(R) ::= NAME(name) LPAREN parlist(list) RPAREN DECLARIZE expression(expr) SEMICOLON . {
  R = ast_define_function(ast, name, list, expr);
}

%type parlist { namelist_t * }
parlist(R) ::= . {
  R = NULL;
}
parlist(R) ::= namelist(list) . {
  R = list;
}

%type namelist { namelist_t * }
namelist(R) ::= NAME(name) . {
  R = ast_new_namelist(ast);
  R = ast_namelist_append(ast, R, name);
}
namelist(R) ::= namelist(sym) COMMA NAME(name) . {
  R = ast_namelist_append(ast, sym, name);
}

%type primary_expression { expr_t * }
primary_expression(R) ::= TRUE|FALSE|NUMBER(tok) . {
  R = ast_new_primary_expr(ast, tok);
}
primary_expression(R) ::= prefix_expression(expr) . {
  R = expr;
}

%type var { expr_t * }
var(R) ::= NAME(name) . {
  R = ast_new_variable(ast, name);
}
//var ::= prefix_expression LBRACKET expression RBRACKET .

%type prefix_expression { expr_t * }
prefix_expression(R) ::= var(var) . {
  R = var;
}
prefix_expression(R) ::= function_call(call) . {
  R = call;
}
prefix_expression(R) ::= OPEN expression(expr) RPAREN . {
  R = ast_new_paren_expr(ast, expr);
}

%type function_call { expr_t * }
function_call(R) ::= prefix_expression(fn) args(args) . {
  R = ast_new_call_expr(ast, fn, args);
}

%type args { expr_t * }
args(R) ::= LPAREN RPAREN . {
  R = NULL;
}
args(R) ::= LPAREN expression_list(list) RPAREN . {
  R = list;
}

%type expression_list { expr_t * }
expression_list(R) ::= expression(expr) . {
  R = expr;
}
expression_list(R) ::= expression_list(list) COMMA expression(expr) . {
  R = ast_expr_comma_append(ast, list, expr);
}

%type unary_expression { expr_t * }
unary_expression(R) ::= primary_expression(expr) . {
  R = expr;
}
unary_expression(R) ::= NOT|MINUS|PLUS(tok) unary_expression(expr) . {
  R = ast_new_unary_expr(ast, tok, expr);
}

%type power_expression { expr_t * }
power_expression(R) ::= unary_expression(expr) . {
  R = expr;
}
power_expression (R) ::= power_expression(left) POWER unary_expression(right) . {
  R = ast_new_pow_expr(ast, left, right);
}

%type multiplicative_expression { expr_t * }
multiplicative_expression(R) ::= power_expression(expr) . {
  R = expr;
}
multiplicative_expression(R) ::= multiplicative_expression(left) TIMES|DIVIDE(tok) power_expression(right) . {
  R = ast_new_mult_expr(ast, left, tok, right);
}

%type additive_expression { expr_t * }
additive_expression(R) ::= multiplicative_expression(expr) . {
  R = expr;
}
additive_expression(R) ::= additive_expression(left) PLUS|MINUS(tok) multiplicative_expression(right) . {
  R = ast_new_add_expr(ast, left, tok, right);
}

%type relational_expression { expr_t * }
relational_expression(R) ::= additive_expression(expr) . {
  R = expr;
}
relational_expression(R) ::= relational_expression(left) LT|LE|GT|GE(tok) additive_expression(right) . {
  R = ast_new_rel_expr(ast, left, tok, right);
}

%type equality_expression { expr_t * }
equality_expression(R)  ::= relational_expression(expr) . {
  R = expr;
}
equality_expression(R)  ::= equality_expression(left) EQ|NE(tok) relational_expression(right) . {
  R = ast_new_eq_expr(ast, left, tok, right);
}

%type logical_and_expression { expr_t * }
logical_and_expression(R) ::= equality_expression(expr) . {
  R = expr;
}
logical_and_expression(R) ::= logical_and_expression(left) AND equality_expression(right) . {
  R = ast_new_and_expr(ast, left, right);
}

%type logical_or_expression { expr_t * }
logical_or_expression(R) ::= logical_and_expression(expr) . {
  R = expr;
}
logical_or_expression(R) ::= logical_or_expression(left) OR logical_and_expression(right) . {
  R = ast_new_or_expr(ast, left, right);
}

%type expression { expr_t * }
expression(R) ::= logical_or_expression(expr) . {
  R = expr;
}
