#ifndef lint
static char rcsid[] = "$Header: /home1/crhet/julian/HYPODD/hypoDD/RCS/datetime_.c,v 1.1 2001/02/15 21:34:55 julian Exp $";
#endif /* lint */

#include <string.h>
#include <time.h>
#include "f77types.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

/*
 * Give current date & time as a character string
 * Full string requires 25 bytes, including null terminator.  Else truncated.
 * f77-callable
 */
void
datetime_(s, n)
	char	*s;
	F77SLA	n;
{
	time_t	now;

	now = time(NULL);
	strncpy(s, ctime(&now), MIN(n-1, 24));
}
