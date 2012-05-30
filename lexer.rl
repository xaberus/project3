#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "expr.h"
#include "expr_parser.h"

%%{
  machine expr;

  count_newline = '\n' @{
    line++;
  };
  any_count_newline = any | count_newline;

  exponent = [eE] [+\-]? digit+;
  number = [1-9] digit* | digit+ '.' digit* exponent? | digit+ exponent | '.' digit+ exponent?;

  alnum_u = alnum | '_';
  alpha_u = alpha | '_';

  main := |*
    ':='  { t = token_new(TOK_DECLARIZE, line, ts, te); fbreak; };

    '('   { t = token_new(TOK_LPAREN, line, ts, te); fbreak; };
    ')'   { t = token_new(TOK_RPAREN, line, ts, te); fbreak; };
    ';'   { t = token_new(TOK_SEMICOLON, line, ts, te); fbreak; };
    ','   { t = token_new(TOK_COMMA, line, ts, te); fbreak; };
    '!'   { t = token_new(TOK_NOT, line, ts, te); fbreak; };
    '='   { t = token_new(TOK_ASSIGN, line, ts, te); fbreak; };

    '^'   { t = token_new(TOK_POWER, line, ts, te); fbreak; };

    '*'   { t = token_new(TOK_TIMES, line, ts, te); fbreak; };
    '/'   { t = token_new(TOK_DIVIDE, line, ts, te); fbreak; };
    'and' { t = token_new(TOK_AND, line, ts, te); fbreak; };

    '+'   { t = token_new(TOK_PLUS, line, ts, te); fbreak; };
    '-'   { t = token_new(TOK_MINUS, line, ts, te); fbreak; };
    'or'  { t = token_new(TOK_OR, line, ts, te); fbreak; };

    '!='  { t = token_new(TOK_NE, line, ts, te); fbreak; };
    '<='  { t = token_new(TOK_LE, line, ts, te); fbreak; };
    '>='  { t = token_new(TOK_GE, line, ts, te); fbreak; };
    '<'   { t = token_new(TOK_LT, line, ts, te); fbreak; };
    '>'   { t = token_new(TOK_GT, line, ts, te); fbreak; };

    any_count_newline - 0x21..0x7E;

    '//' [^\n]* count_newline;
    '/*' any_count_newline* :>> '*/';

    number { t = token_new(TOK_NUMBER, line, ts, te); fbreak; };

    alpha_u alnum_u* { t = token_new(TOK_NAME, line, ts, te); fbreak;};
  *|;
}%%

token_t * token_new(int type, int line, const char * start, const char * end) {
  token_t * t = malloc(sizeof(token_t)); assert(t);
  t->type = type;
  t->line = line;
  t->start = start;
  t->end = end;
  t->next = NULL;
  return t;
}

token_t * exprlexer_tokenize(const char * buf) {
  %% write data;
  token_t * tokens = NULL, * t = NULL;
  token_t ** append = &tokens;
  int cs, line = 1, act;
  const char * ts, * te;
  const char * p = buf, * pe = buf + strlen(buf);
  const char * eof = pe;
  %% write init;
  while (p < pe) {
    %% write exec;
    if (cs == expr_error) {
      fprintf(stderr, "syntax error at %s\n", p);
      break;
    }
    if (t) {
      *append = t;
      append = &t->next;
    }
  }
  return tokens;
}