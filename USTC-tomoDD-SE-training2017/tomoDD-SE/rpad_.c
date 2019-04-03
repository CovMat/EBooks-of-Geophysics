#ifndef lint
static char rcsid[] = "$Header: /home1/crhet/julian/HYPODD/hypoDD/RCS/rpad_.c,v 1.1 2001/02/15 20:02:38 julian Exp $";
#endif /* lint */

#include <ctype.h>
#include "f77types.h"

/*
 * Pad character string to full length with white space (remove null terminator)
 * f77-callable
 */
void
rpad_(s, n)
	char	*s;
	F77SLA	n;
{
	int	i, j;

	for (i=0; i<n; i++) {
		if (s[i] == '\0') {
			/* Found terminator: Put spaces from here to end */
			for (j=i; j<n; j++)
				s[j] = ' ';
			break;
		}
	}
}
