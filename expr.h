#ifndef _EXPR_H_
#define _EXPR_H_

#include <stdint.h>
#include <complex.h>
#include "hashmap.h"
/*

  L = 10;
  V = 1000;

  potential(x, t) := if abs(x) < L/2 then V else 0;
  psi(x) := exp(-x*x/(2)) * cexp(I * 4 * x);

*/

/*typedef enum {
  VALUE_FLAG_TRUE       = 1 << 0,
  VALUE_FLAG_FALSE      = 1 << 1,
  VALUE_FLAG_BOOL       = VALUE_FLAG_TRUE | VALUE_FLAG_FALSE,
  VALUE_FLAG_NUMBER     = 1 << 2,
} expr_flag_t;
*/

typedef struct {
  complex double number;
} value_t;

typedef enum {
  EXPR_NONE,

  EXPR_TRUE,
  EXPR_FALSE,

  EXPR_NUMBER,
  EXPR_VARIABLE,

  EXPR_PAREN,

  EXPR_POW,

  EXPR_ADD,
  EXPR_SUB,
  EXPR_MULT,
  EXPR_DIV,

  EXPR_LT,
  EXPR_LE,
  EXPR_GT,
  EXPR_GE,
  EXPR_EQ,
  EXPR_NE,

  EXPR_NEG,
  EXPR_NOT,
  EXPR_AND,
  EXPR_OR,

  EXPR_FN,
} expr_type_t;

typedef struct token token_t;
typedef struct ast ast_t;
typedef struct deflist deflist_t;
typedef struct symdef symdef_t;
typedef struct vardef vardef_t;
typedef struct name name_t;
typedef struct namelist namelist_t;
typedef struct funcdef funcdef_t;
typedef struct expr expr_t;

struct expr {
  unsigned  type;
  value_t * value;
  name_t  * leash; /* btw., where is the hound? */
  char    * name;
  expr_t  * left;
  expr_t  * right;
  expr_t  * comma;

  expr_t  * next; // global list
};


struct token {
  int          type;
  int          line;
  const char * start;
  const char * end;
  token_t    * next;
};

#define EXPR_SYM_VAR  1
#define EXPR_SYM_FUNC 2

struct symdef {
  int        type;
  char     * name;
  symdef_t * next;
};


struct vardef {
  symdef_t   sym;
  expr_t   * value;
};

struct funcdef {
  symdef_t     sym;
  namelist_t * params;
  expr_t     * value;
};

struct name {
  char    * name;
  value_t * value;
};

struct namelist {
  int      length;
  int      alloc;
  name_t   names[];
};

struct deflist {
  symdef_t * list;
  symdef_t ** append;
};

typedef struct {
  int    alloc;
  int    length;
  char * data;
} pbuf_t;

struct ast {
  token_t    * tokens;
  hashmap_t  * symbols;
  symdef_t   * sdefs;
  symdef_t  ** sappend;
  expr_t     * exprs;
  hashmap_t  * values;
  pbuf_t     * buf;
};

ast_t * ast_new(token_t * tokens);

void ast_free(ast_t * ast);

token_t * exprlexer_tokenize(const char * buf);

void ast_delete_token(ast_t * ast, token_t * tok);

expr_t * ast_expr_clone(ast_t * ast, expr_t * orig);

vardef_t * ast_define_variable(ast_t * ast, token_t * name, expr_t * expr);
void ast_deflist_append_var(ast_t * ast, vardef_t * var);
namelist_t * ast_new_namelist(ast_t * ast);
namelist_t * ast_namelist_append(ast_t * ast, namelist_t * list, token_t * name);
funcdef_t * ast_define_function(ast_t * ast, token_t * name, namelist_t * params, expr_t * expr);
void ast_deflist_append_func(ast_t * ast, funcdef_t * func);
expr_t * ast_expr_comma_append(ast_t * ast, expr_t * list, expr_t * expr);
expr_t * ast_new_primary_expr(ast_t * ast, token_t * tok);
expr_t * ast_new_variable(ast_t * ast, token_t * name);
expr_t * ast_new_paren_expr(ast_t * ast, expr_t * expr);
expr_t * ast_new_call_expr(ast_t * ast, expr_t * fn, expr_t * args);
expr_t * ast_new_unary_expr(ast_t * ast, token_t * tok, expr_t * expr);
expr_t * ast_new_pow_expr(ast_t * ast, expr_t * left, expr_t * right);
expr_t * ast_new_mult_expr(ast_t * ast, expr_t * left, token_t * tok, expr_t * right);
expr_t * ast_new_add_expr(ast_t * ast, expr_t * left, token_t * tok, expr_t * right);
expr_t * ast_new_rel_expr(ast_t * ast, expr_t * left, token_t * tok, expr_t * right);
expr_t * ast_new_eq_expr(ast_t * ast, expr_t * left, token_t * tok, expr_t * right);
expr_t * ast_new_and_expr(ast_t * ast, expr_t * left, expr_t * right);
expr_t * ast_new_or_expr(ast_t * ast, expr_t * left, expr_t * right);

#endif /* _EXPR_H_ */