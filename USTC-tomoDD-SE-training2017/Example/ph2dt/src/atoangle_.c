#ifndef lint
static char rcsid[]="$Header: /home1/crhet/julian/HYPODD/ph2dt/RCS/atoangle_.c,v 1.1 2001/02/13 21:03:41 julian Exp $";
#endif lint
#include <string.h>
#include "compat.h"
#include "f77types.h"

	double	atoangle PARMS((char*));

/*
 * f77-callable interface to atoangle function
 */
double
atoangle_(p, l)
	char	*p;	/* Character string		(input)	*/
	F77SLA	l;	/* Length of *p			(input)	*/
{
	char	q[l+1];	/* Buffer, with space for terminator	*/

	return atoangle(strncpy(q, p, l));
}
