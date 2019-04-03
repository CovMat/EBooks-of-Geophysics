#ifndef lint
static char rcsid[] = "$Header: /home1/crhet/julian/HYPODD/hypoDD/RCS/hypot_.c,v 1.1 2001/02/14 00:21:19 julian Exp $";
#endif /* lint */

#include <math.h>

/* f77-callable interface to hypot function */
double
hypot_(a, b)
	float	*a, *b;
{
	return hypot((double)*a, (double)*b);
}
