#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include "expr.h"
#include "expr_parser.h"

void *expr_parserAlloc(void *(*mallocProc)(size_t));
void *expr_parserAlloc(void *(*mallocProc)(size_t));
void expr_parserFree(void * p, void (*freeProc)(void*));

void expr_parser(void *p, int type, token_t * tok, ast_t * ast);

/***********************************************/

pbuf_t * pbuf_new(int alloc) {
  pbuf_t * b = malloc(sizeof(pbuf_t)); assert(b);
  b->data = malloc(alloc); assert(b->data);
  b->alloc = alloc;
  b->length = 0;
  return b;
}

void pbuf_reset(pbuf_t * b) {
  b->length = 0;
}

void pbuf_free(pbuf_t * b) {
  free(b->data);
  free(b);
}


void pbuf_append(pbuf_t * b, char c) {
  if (b->length < b->alloc) {
    b->data[b->length++] = c;
  } else {
    int alloc = b->alloc * 2;
    char * t = realloc(b->data, alloc); assert(t);
    b->data = t;
    b->data[b->length++] = c;
  }
}

void pbuf_appends(pbuf_t * b, char * s) {
  int len = strlen(s);
  if (b->length + len < b->alloc) {
    for (; *s; s++) {
      b->data[b->length++] = *s;
    }
  } else {
    int alloc = (b->alloc + len) * 2;
    char * t = realloc(b->data, alloc); assert(t);
    b->data = t;
    for (; *s; s++) {
      b->data[b->length++] = *s;
    }
  }
}

void pbuf_appendf(pbuf_t * b, const char * fmt, ...) {
  char buf[512];

  va_list ap;

  va_start(ap, fmt);
  vsnprintf(buf, sizeof(buf), fmt, ap);
  va_end(ap);

  pbuf_appends(b, buf);
}

/***********************************************/

ast_t * ast_new(token_t * tokens) {
  token_t * t = tokens;
  ast_t * ast = malloc(sizeof(ast_t)); assert(ast);

  //expr_parserTrace(stderr, ">>>> ");

  ast->tokens = tokens;

  ast->symbols = hashmap_new(255);
  ast->values = hashmap_new(255);

  ast->exprs = NULL;
  ast->sdefs = NULL;
  ast->sappend = &ast->sdefs;

  void * parser = expr_parserAlloc(malloc); assert(parser);
  while (t) {
    //printf("TOK %d <%s>\n", t->type, t->start ? t->value : "");
    expr_parser(parser, t->type, t, ast);
    t = t->next;
  }
  expr_parser(parser, 0, NULL, ast);
  expr_parserFree(parser, free);

  return ast;
}

funcdef_t * ast_get_func(ast_t * ast, const char * fn) {
  symdef_t * s = hashmap_get(ast->symbols, fn);
  if (!s) {
    fprintf(stderr, "function '%s' ist not defined\n", fn);
    abort();
  }
  if (s->type != EXPR_SYM_FUNC) {
    fprintf(stderr, "symbol '%s' ist not a function\n", fn);
    abort();
  }
  return (funcdef_t *) s;
}

complex double ast_eval(ast_t * ast, funcdef_t fdef, int argc, complex double argv) {
  symdef_t * s = hashmap_get(ast->symbols, fn);
  if (!s) {
    fprintf(stderr, "function '%s' ist not defined\n", fn);
    abort();
  }
  if (s->type != EXPR_SYM_FUNC) {
    fprintf(stderr, "symbol '%s' ist not a function\n", fn);
    abort();
  }
  return (funcdef_t *) s;
}


void expr_free(expr_t * expr) {
  switch (expr->type) {
    case EXPR_NONE:
    case EXPR_TRUE:
    case EXPR_FALSE:
    case EXPR_NUMBER:
    case EXPR_PAREN:
    case EXPR_POW:
    case EXPR_ADD:
    case EXPR_SUB:
    case EXPR_MULT:
    case EXPR_DIV:
    case EXPR_LT:
    case EXPR_LE:
    case EXPR_GT:
    case EXPR_GE:
    case EXPR_EQ:
    case EXPR_NE:
    case EXPR_NEG:
    case EXPR_NOT:
    case EXPR_AND:
    case EXPR_OR:
      break;
    case EXPR_VARIABLE:
    case EXPR_FN:
      free(expr->name);
      break;
  }
  free(expr);
}

void * value_free(const char * key, void * o) {
  free(o);
  return NULL;
}

void ast_free(ast_t * ast) {
  hashmap_free(ast->symbols);

  hashmap_foreach(ast->values, value_free);
  hashmap_free(ast->values);

  expr_t * e = ast->exprs;
  while (e) {
    expr_t * n = e->next;
    expr_free(e);
    e = n;
  }

  token_t * t = ast->tokens;
  while (t) {
    token_t * n = t->next;
    free(t);
    t = n;
  }

  //assert(ast->defs);
  symdef_t * s = ast->sdefs;
  while (s) {
    symdef_t * n = s->next;
    free(s->name);
    switch (s->type) {
      case EXPR_SYM_VAR: {
        //vardef_t * vdef = (vardef_t *) s;
      } break;
      case EXPR_SYM_FUNC: {
        funcdef_t * fdef = (funcdef_t *) s;
        namelist_t * params = fdef->params;
        for (int k = 0; k < params->length; k++) {
          free(params->names[k].name);
        }
        free(fdef->params);
      } break;
    }
    free(s);
    s = n;
  }

  free(ast);
}

vardef_t * ast_define_variable(ast_t * ast, token_t * name, expr_t * expr) {
  assert(name->type == TOK_NAME);
  assert(name->start);
  assert(expr);

  char * key = strndup(name->start, name->end - name->start); assert(key);

  symdef_t * s = hashmap_get(ast->symbols, key);
  if (s) {
    fprintf(stderr, "warning: redefining symbol '%s'\n", key);
  }

  vardef_t * vdef = malloc(sizeof(vardef_t)); assert(vdef);
  vdef->sym.type = EXPR_SYM_VAR;
  vdef->sym.name = key;
  vdef->sym.next = NULL;
  vdef->value = expr;

  return vdef;
}

void ast_deflist_append_var(ast_t * ast, vardef_t * var) {
  *ast->sappend = &var->sym;
  ast->sappend = &var->sym.next;
  hashmap_insert(ast->symbols, var->sym.name, var);
}

namelist_t * ast_new_namelist(ast_t * ast) {
  (void) ast;
  namelist_t * names = malloc(sizeof(namelist_t) + sizeof(name_t) * 4); assert(names);
  names->alloc = 4;
  names->length = 0;
  return names;
}

namelist_t * ast_namelist_append(ast_t * ast, namelist_t * list, token_t * name) {
  (void) ast;

  char * key = strndup(name->start, name->end - name->start); assert(key);

  if (list->length < list->alloc) {
    list->names[list->length].name = key;
    list->names[list->length++].value = NULL;
    return list;
  } else {
    int alloc = list->alloc * 2;
    namelist_t * t = realloc(list, sizeof(namelist_t) + sizeof(name_t) * alloc); assert(t);
    t->alloc = alloc;
    t->names[t->length].name = key;
    t->names[t->length++].value = NULL;
    return t;
  }
}

/***********************************************/

void expr_tostring(expr_t * expr, pbuf_t * buf);

void expr_tostring_i(expr_t * expr, pbuf_t * buf) {
  if (expr->value) {
    complex double z = expr->value->number;
    double x = creal(z), y = cimag(z);
    if (y != 0) {
      pbuf_appendf(buf, "(%g + I %g)", x, y);
    } else {
      pbuf_appendf(buf, "%g", x);
    }
    return;
  }
  switch (expr->type) {
    case EXPR_NONE:
      pbuf_appends(buf, "<none>");
      break;
    case EXPR_TRUE:
      pbuf_appends(buf, "true");
      break;
    case EXPR_FALSE:
      pbuf_appends(buf, "false");
      break;
    case EXPR_NUMBER: {
      complex double z = expr->value->number;
      double x = creal(z), y = cimag(z);
      if (y != 0) {
        pbuf_appendf(buf, "%g + I %g", x, y);
      } else {
        pbuf_appendf(buf, "%g", x);
      }
    } break;
    case EXPR_PAREN:
      pbuf_appends(buf, "(");
      expr_tostring(expr->left, buf);
      pbuf_appends(buf, ")");
      break;
    case EXPR_POW:
      expr_tostring(expr->left, buf);
      pbuf_appends(buf, "^");
      expr_tostring(expr->right, buf);
      break;
    case EXPR_ADD:
      expr_tostring(expr->left, buf);
      pbuf_appends(buf, "+");
      expr_tostring(expr->right, buf);
      break;
    case EXPR_SUB:
      expr_tostring(expr->left, buf);
      pbuf_appends(buf, "-");
      expr_tostring(expr->right, buf);
      break;
    case EXPR_MULT:
      expr_tostring(expr->left, buf);
      pbuf_appends(buf, "*");
      expr_tostring(expr->right, buf);
      break;
    case EXPR_DIV:
      expr_tostring(expr->left, buf);
      pbuf_appends(buf, "/");
      expr_tostring(expr->right, buf);
      break;
    case EXPR_LT:
      expr_tostring(expr->left, buf);
      pbuf_appends(buf, "<");
      expr_tostring(expr->right, buf);
      break;
    case EXPR_LE:
      expr_tostring(expr->left, buf);
      pbuf_appends(buf, "<=");
      expr_tostring(expr->right, buf);
      break;
    case EXPR_GT:
      expr_tostring(expr->left, buf);
      pbuf_appends(buf, ">");
      expr_tostring(expr->right, buf);
      break;
    case EXPR_GE:
      expr_tostring(expr->left, buf);
      pbuf_appends(buf, ">=");
      expr_tostring(expr->right, buf);
      break;
    case EXPR_EQ:
      expr_tostring(expr->left, buf);
      pbuf_appends(buf, "==");
      expr_tostring(expr->right, buf);
      break;
    case EXPR_NE:
      expr_tostring(expr->left, buf);
      pbuf_appends(buf, "!=");
      expr_tostring(expr->right, buf);
      break;
    case EXPR_NEG:
      pbuf_appends(buf, "-");
      expr_tostring(expr->left, buf);
      break;
    case EXPR_NOT:
      pbuf_appends(buf, "not");
      expr_tostring(expr->left, buf);
      break;
    case EXPR_AND:
      pbuf_appends(buf, "and");
      expr_tostring(expr->left, buf);
      break;
    case EXPR_OR:
      pbuf_appends(buf, "or");
      expr_tostring(expr->left, buf);
      break;
    case EXPR_VARIABLE:
      if (expr->leash) {
        pbuf_appendf(buf, "[%s]", expr->leash->name);
      } else {
        pbuf_appendf(buf, "%s", expr->name);
      }
      break;
    case EXPR_FN: {
      expr_tostring(expr->left, buf);
      pbuf_appends(buf, "(");
      if (expr->right) {
        expr_tostring(expr->right, buf);
      }
      pbuf_appends(buf, ")");
    } break;
  }
}

void expr_tostring(expr_t * expr, pbuf_t * buf) {
  if (expr->comma) {
    pbuf_append(buf, '(');
    expr_t * b = expr;
    while (b) {
      expr_tostring_i(b, buf);
      b = b->comma;
      if (b) {
        pbuf_appends(buf, ", ");
      }
    }
    pbuf_append(buf, ')');
  } else {
    expr_tostring_i(expr, buf);
  }
}

void expr_print(expr_t * expr) {
  pbuf_t * b = pbuf_new(512);
  expr_tostring(expr, b);
  printf("%.*s\n", b->length, b->data);
  pbuf_free(b);
}

/***********************************************/

expr_t * expr_mark_vars(ast_t * ast, namelist_t * local, expr_t * expr) {
  expr_t ** p = &expr, * b = expr;
  while (b) {
    switch (b->type) {
      case EXPR_NONE:
      case EXPR_TRUE:
      case EXPR_FALSE:
      case EXPR_NUMBER:
        break;
      case EXPR_PAREN:
        b->left = expr_mark_vars(ast, local, b->left);
        break;
      case EXPR_POW:
      case EXPR_ADD:
      case EXPR_SUB:
      case EXPR_MULT:
      case EXPR_DIV:
      case EXPR_LT:
      case EXPR_LE:
      case EXPR_GT:
      case EXPR_GE:
      case EXPR_EQ:
      case EXPR_NE:
      case EXPR_AND:
      case EXPR_OR:
        b->left = expr_mark_vars(ast, local, b->left);
        b->right = expr_mark_vars(ast, local, b->right);
        break;
      case EXPR_NEG:
      case EXPR_NOT:
        b->left = expr_mark_vars(ast, local, b->left);
        break;
      case EXPR_VARIABLE: {
        int hound = -1;
        for (int k = 0; k < local->length; k++) {
          if (strcmp(local->names[k].name, b->name) == 0) {
            hound = k;
            break;
          }
        }
        if (hound >= 0) {
          b->leash = &local->names[hound];
        } else {
          /* global variable found? */
          symdef_t * s = hashmap_get(ast->symbols, b->name);
          if (!s) {
            fprintf(stderr, "error: undefined symbol '%s'\n", b->name);
            abort();
          } else {
            switch (s->type) {
              case EXPR_SYM_VAR: {
                expr_t * n = ast_expr_clone(ast, ((vardef_t *) s)->value);
                n = ast_new_paren_expr(ast, n);
                *p = n;
                n->comma = b->comma;
                b = n;
                continue;
              } break;
              case EXPR_SYM_FUNC: {
                fprintf(stderr, "error: cannot use '%s' as variable\n", b->name);
                abort();
              } break;
            }
          }
        }
      } break;
      case EXPR_FN:
        b->right = expr_mark_vars(ast, local, b->right);
        break;
    }
    p = &b->comma;
    b = *p;
  }
  return expr;
}

/**************************************************************/

value_t * value_pow(value_t * a, value_t * b) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = cpow(a->number, b->number);
  return val;
}

value_t * value_add(value_t * a, value_t * b) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = a->number + b->number;
  return val;
}

value_t * value_sub(value_t * a, value_t * b) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = a->number - b->number;
  return val;
}

value_t * value_mult(value_t * a, value_t * b) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = a->number * b->number;
  return val;
}

value_t * value_div(value_t * a, value_t * b) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = a->number / b->number;
  return val;
}

value_t * value_lt(value_t * a, value_t * b) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = creal(a->number) < creal(b->number);
  return val;
}

value_t * value_le(value_t * a, value_t * b) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = creal(a->number) <= creal(b->number);
  return val;
}

value_t * value_gt(value_t * a, value_t * b) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = creal(a->number) > creal(b->number);
  return val;
}

value_t * value_ge(value_t * a, value_t * b) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = creal(a->number) >= creal(b->number);
  return val;
}

value_t * value_eq(value_t * a, value_t * b) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = a->number == b->number;
  return val;
}

value_t * value_ne(value_t * a, value_t * b) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = a->number != b->number;
  return val;
}

value_t * value_and(value_t * a, value_t * b) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = (a->number == 1) && (b->number == 1);
  return val;
}

value_t * value_or(value_t * a, value_t * b) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = (a->number == 1) || (b->number == 1);
  return val;
}

value_t * value_neg(value_t * a) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = -a->number;
  return val;
}

value_t * value_not(value_t * a) {
  value_t * val = malloc(sizeof(expr_t)); assert(val);
  val->number = !(a->number == 1);
  return val;
}

/**************************************************************/

void expr_fold(ast_t * ast, expr_t * expr, pbuf_t * buf);

int expr_release(expr_t * expr) {
  expr_t * b = expr;
  while (b) {
    switch (b->type) {
      case EXPR_NONE:
      case EXPR_TRUE:
      case EXPR_FALSE:
      case EXPR_NUMBER:
        break;
      case EXPR_NEG:
      case EXPR_NOT:
      case EXPR_PAREN:
        if (expr_release(b->left)) { b->value = NULL; return 1; }
        break;
      case EXPR_POW:
      case EXPR_ADD:
      case EXPR_SUB:
      case EXPR_MULT:
      case EXPR_DIV:
      case EXPR_LT:
      case EXPR_LE:
      case EXPR_GT:
      case EXPR_GE:
      case EXPR_EQ:
      case EXPR_NE:
      case EXPR_AND:
      case EXPR_OR:
        if (expr_release(b->left)) { b->value = NULL; return 1; }
        if (expr_release(b->right)) { b->value = NULL; return 1; }
        break;
      case EXPR_VARIABLE:
        if (b->leash) {
          b->value = NULL;
          return 1;
        }
      case EXPR_FN:
        if (expr_release(b->right)) { b->value = NULL; return 1; }
        break;
    }
    b = b->comma;
  }
  return 0;
}

value_t * expr_call_fold(ast_t * ast, expr_t * expr, pbuf_t * buf) {
  assert(expr->type = EXPR_FN);
  expr_t * fn = expr->left;
  assert(fn->type == EXPR_VARIABLE);
  symdef_t * s = hashmap_get(ast->symbols, expr->left->name);
  assert(s && s->type == EXPR_SYM_FUNC);
  funcdef_t * fdef = (funcdef_t *) s;
  expr_t * args = expr->right;

  namelist_t * params = fdef->params;

  int k = 0;
  while (args) {
    assert(k < params->length);
    params->names[k++].value = args->value;
    args = args->comma;
  }

  expr_fold(ast, fdef->value, buf);

  assert(fdef->value->value);

  expr->value = fdef->value->value;

  for (k = 0; k < params->length; k++) {
    params->names[k++].value = NULL;
  }

  value_t * val = fdef->value->value;
  assert(val);

  if (expr_release(fdef->value)) {
    fdef->value->value = NULL;
  }

  return val;
}

void expr_fold(ast_t * ast, expr_t * expr, pbuf_t * buf) {

#define fold_eval(_expr, _val)\
  do {\
    pbuf_reset(buf);\
    expr_tostring((_expr), buf);\
    pbuf_append(buf, '\0');\
    value_t * val = hashmap_get(ast->values, buf->data);\
    if (!val) {\
      val = (_val);\
      hashmap_insert(ast->values, buf->data, val);\
    }\
    (_expr)->value = val;\
  } while (0)

  expr_t * b = expr;
  while (b) {
    switch (b->type) {
      case EXPR_NONE:
      case EXPR_TRUE:
      case EXPR_FALSE:
      case EXPR_NUMBER:
        break;
      case EXPR_PAREN:
        expr_fold(ast, b->left, buf);
        if (b->left->value) {
          b->value = b->left->value;
        }
        break;
      case EXPR_POW:
        expr_fold(ast, b->left, buf); expr_fold(ast, b->right, buf);
        if (b->left->value && b->right->value) {
          fold_eval(b, value_pow(b->left->value, b->right->value));
        }
        break;
      case EXPR_ADD:
        expr_fold(ast, b->left, buf); expr_fold(ast, b->right, buf);
        if (b->left->value && b->right->value) {
          fold_eval(b, value_add(b->left->value, b->right->value));
        }
        break;
      case EXPR_SUB:
        expr_fold(ast, b->left, buf); expr_fold(ast, b->right, buf);
        if (b->left->value && b->right->value) {
          fold_eval(b, value_sub(b->left->value, b->right->value));
        }
        break;
      case EXPR_MULT:
        expr_fold(ast, b->left, buf); expr_fold(ast, b->right, buf);
        if (b->left->value && b->right->value) {
          fold_eval(b, value_mult(b->left->value, b->right->value));
        }
        break;
      case EXPR_DIV:
        expr_fold(ast, b->left, buf); expr_fold(ast, b->right, buf);
        if (b->left->value && b->right->value) {
          fold_eval(b, value_div(b->left->value, b->right->value));
        }
        break;
      case EXPR_LT:
        expr_fold(ast, b->left, buf); expr_fold(ast, b->right, buf);
        if (b->left->value && b->right->value) {
          fold_eval(b, value_lt(b->left->value, b->right->value));
        }
        break;
      case EXPR_LE:
        expr_fold(ast, b->left, buf); expr_fold(ast, b->right, buf);
        if (b->left->value && b->right->value) {
          fold_eval(b, value_le(b->left->value, b->right->value));
        }
        break;
      case EXPR_GT:
        expr_fold(ast, b->left, buf); expr_fold(ast, b->right, buf);
        if (b->left->value && b->right->value) {
          fold_eval(b, value_gt(b->left->value, b->right->value));
        }
        break;
      case EXPR_GE:
        expr_fold(ast, b->left, buf); expr_fold(ast, b->right, buf);
        if (b->left->value && b->right->value) {
          fold_eval(b, value_ge(b->left->value, b->right->value));
        }
        break;
      case EXPR_EQ:
        expr_fold(ast, b->left, buf); expr_fold(ast, b->right, buf);
        if (b->left->value && b->right->value) {
          fold_eval(b, value_eq(b->left->value, b->right->value));
        }
        break;
      case EXPR_NE:
        expr_fold(ast, b->left, buf); expr_fold(ast, b->right, buf);
        if (b->left->value && b->right->value) {
          fold_eval(b, value_ne(b->left->value, b->right->value));
        }
        break;
      case EXPR_AND:
        expr_fold(ast, b->left, buf); expr_fold(ast, b->right, buf);
        if (b->left->value && b->right->value) {
          fold_eval(b, value_and(b->left->value, b->right->value));
        }
        break;
      case EXPR_OR:
        expr_fold(ast, b->left, buf); expr_fold(ast, b->right, buf);
        if (b->left->value && b->right->value) {
          fold_eval(b, value_or(b->left->value, b->right->value));
        }
        break;
      case EXPR_NEG:
        expr_fold(ast, b->left, buf);
        if (b->left->value) {
          fold_eval(b, value_neg(b->left->value));
        }
        break;
      case EXPR_NOT:
        expr_fold(ast, b->left, buf);
        if (b->left->value) {
          fold_eval(b, value_not(b->left->value));
        }
        break;
      case EXPR_VARIABLE:
        if (b->leash && b->leash->value) {
          b->value = b->leash->value;
        }
        break;
      case EXPR_FN: {
        expr_fold(ast, b->right, buf);
        int eval = 1;
        expr_t * t = b->right;
        while (t) {
          if (!t->value) {
            eval = 0;
          }
          t = t->comma;
        }
        if (eval) {
          fold_eval(b, expr_call_fold(ast, b, buf));
        }
      } break;
    }
    b = b->comma;
  }
}

expr_t * ast_fold_expr(ast_t * ast, namelist_t * local, expr_t * expr) {
  expr = expr_mark_vars(ast, local, expr);
  pbuf_t * b = pbuf_new(512);
  expr_fold(ast, expr, b);
  pbuf_free(b);
  expr_print(expr);
  return expr;
}

funcdef_t * ast_define_function(ast_t * ast, token_t * name, namelist_t * params, expr_t * expr) {
  assert(name->type == TOK_NAME);
  assert(name->start);
  assert(expr);

  expr = ast_fold_expr(ast, params, expr);

  char * key = strndup(name->start, name->end - name->start); assert(key);

  symdef_t * s = hashmap_get(ast->symbols, key);
  if (s) {
    fprintf(stderr, "warning: redefining symbol '%s'\n", key);
  }

  funcdef_t * fdef = malloc(sizeof(funcdef_t)); assert(fdef);
  fdef->sym.type = EXPR_SYM_FUNC;
  fdef->sym.name = key;
  fdef->sym.next = NULL;
  fdef->params = params;
  fdef->value = expr;

  return fdef;
}

void ast_deflist_append_func(ast_t * ast, funcdef_t * func) {
  *ast->sappend = &func->sym;
  ast->sappend = &func->sym.next;
  hashmap_insert(ast->symbols, func->sym.name, func);
}


expr_t * ast_expr_comma_append(ast_t * ast, expr_t * list, expr_t * expr) {
  (void) ast;

  expr_t ** p = &list->comma, * b = *p;
  while (b) {
    p = &b->comma;
    b = *p;
  }
  *p = expr;

  return list;
}

expr_t * ast_new_empty_expr(ast_t * ast) {
  expr_t * expr = malloc(sizeof(expr_t)); assert(expr);
  expr->type = EXPR_NONE;
  expr->name = NULL;
  expr->comma = NULL;
  expr->leash = NULL;
  expr->left = NULL;
  expr->right = NULL;
  expr->next = ast->exprs;
  ast->exprs = expr;
  return expr;
}

expr_t * ast_expr_clone(ast_t * ast, expr_t * orig) {
  expr_t * expr = malloc(sizeof(expr_t)); assert(expr);
  expr->type = orig->type;
  if (orig->name) {
    expr->name = strdup(orig->name); assert(expr->name);
  }
  expr->left = orig->left;
  expr->right = orig->right;
  expr->leash = orig->leash;
  expr->value = orig->value;

  expr->comma = NULL;

  expr->next = ast->exprs;
  ast->exprs = expr;

  return expr;
}

value_t * ast_new_number_value_from_token(ast_t * ast, token_t * tok) {
  double d = strtod(tok->start, NULL);
  char buf[50];
  snprintf(buf, sizeof(buf), "%.30e", d);

  value_t * val = hashmap_get(ast->values, buf);
  if (!val) {
    val = malloc(sizeof(expr_t)); assert(val);
    val->number = d;
    hashmap_insert(ast->values, buf, val);
  }

  return val;
}

expr_t * ast_new_primary_expr(ast_t * ast, token_t * tok) {
  expr_t * expr = ast_new_empty_expr(ast);
  switch (tok->type) {
    case TOK_FALSE:
      expr->type = EXPR_FALSE;
      expr->value = NULL;
      break;
    case TOK_TRUE:
      expr->type = EXPR_TRUE;
      expr->value = NULL;
      break;
    case TOK_NUMBER:
      expr->type = EXPR_NUMBER;
      expr->value = ast_new_number_value_from_token(ast, tok);
      break;
  }
  return expr;
}

expr_t * ast_new_variable(ast_t * ast, token_t * name) {
  expr_t * expr = ast_new_empty_expr(ast);
  char * key = strndup(name->start, name->end - name->start); assert(key);

  expr->type = EXPR_VARIABLE;
  expr->name = key;
  expr->value = NULL;

  return expr;
}

expr_t * ast_new_paren_expr(ast_t * ast, expr_t * paren) {
  expr_t * expr = ast_new_empty_expr(ast);

  expr->type = EXPR_PAREN;
  expr->name = NULL;
  expr->left = paren;
  expr->value = NULL;

  return expr;
}

expr_t * ast_new_call_expr(ast_t * ast, expr_t * fn, expr_t * args) {
  expr_t * expr = ast_new_empty_expr(ast);

  expr->type = EXPR_FN;
  expr->name = NULL;
  expr->left = fn;
  expr->right = args;
  expr->value = NULL;

  return expr;
}

expr_t * ast_new_unary_expr(ast_t * ast, token_t * tok, expr_t * un) {
  expr_t * expr;
  switch (tok->type) {
    case TOK_NOT: {
      expr_t * expr = ast_new_empty_expr(ast);
      expr->type = EXPR_NOT;
      expr->left = un;
      expr->value = NULL;
    } break;
    case TOK_MINUS: {
      expr_t * expr = ast_new_empty_expr(ast);
      expr->type = EXPR_NEG;
      expr->left = un;
      expr->value = NULL;
    } break;
    case TOK_PLUS:
      expr = un;
  }
  return expr;
}

expr_t * ast_new_pow_expr(ast_t * ast, expr_t * left, expr_t * right) {
  expr_t * expr = ast_new_empty_expr(ast);
  expr->type = EXPR_POW;
  expr->left = left;
  expr->right = right;
  expr->value = NULL;
  return expr;
}

expr_t * ast_new_mult_expr(ast_t * ast, expr_t * left, token_t * tok, expr_t * right) {
  expr_t * expr = ast_new_empty_expr(ast);
  switch (tok->type) {
    case TOK_TIMES:
      expr->type = EXPR_MULT;
      break;
    case TOK_DIVIDE:
      expr->type = EXPR_DIV;
      break;
  }
  expr->left = left;
  expr->right = right;
  expr->value = NULL;
  return expr;
}

expr_t * ast_new_add_expr(ast_t * ast, expr_t * left, token_t * tok, expr_t * right) {
  expr_t * expr = ast_new_empty_expr(ast);
  switch (tok->type) {
    case TOK_PLUS:
      expr->type = EXPR_ADD;
      break;
    case TOK_MINUS:
      expr->type = EXPR_SUB;
      break;
  }
  expr->left = left;
  expr->right = right;
  expr->value = NULL;
  return expr;
}
expr_t * ast_new_rel_expr(ast_t * ast, expr_t * left, token_t * tok, expr_t * right) {
  expr_t * expr = ast_new_empty_expr(ast);
  switch (tok->type) {
    case TOK_LT:
      expr->type = EXPR_LT;
      break;
    case TOK_LE:
      expr->type = EXPR_LE;
      break;
    case TOK_GT:
      expr->type = EXPR_GT;
      break;
    case TOK_GE:
      expr->type = EXPR_GE;
      break;
  }
  expr->left = left;
  expr->right = right;
  expr->value = NULL;
  return expr;
}

expr_t * ast_new_eq_expr(ast_t * ast, expr_t * left, token_t * tok, expr_t * right) {
  expr_t * expr = ast_new_empty_expr(ast);
  switch (tok->type) {
    case TOK_EQ:
      expr->type = EXPR_EQ;
      break;
    case TOK_NE:
      expr->type = EXPR_NE;
      break;
  }
  expr->left = left;
  expr->right = right;
  expr->value = NULL;
  return expr;
}
expr_t * ast_new_and_expr(ast_t * ast, expr_t * left, expr_t * right) {
  expr_t * expr = ast_new_empty_expr(ast);
  expr->type = EXPR_AND;
  expr->left = left;
  expr->right = right;
  expr->value = NULL;
  return expr;
}
expr_t * ast_new_or_expr(ast_t * ast, expr_t * left, expr_t * right) {
  expr_t * expr = ast_new_empty_expr(ast);
  expr->type = EXPR_OR;
  expr->left = left;
  expr->right = right;
  expr->value = NULL;
  return expr;
}