
#if !defined(CODERERROR_H)
#define CODERERROR_H

extern int coder_error_flag;

int return_coder_error_flag();

void set_coder_error_flag(int flag);

#define CODER_ERR_NAN        1   /* NaN generated */
#define CODER_ERR_BOUNDS     2   /* Improper bisection bounds */
#define CODER_ERR_MAXITS     3   /* Maximum iterations reached */

#endif
