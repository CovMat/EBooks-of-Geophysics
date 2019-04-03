#ifndef lint
static char rcsid[]="$Header: /home1/crhet/julian/HYPODD/hypoDD/RCS/atoangle_.c,v 1.1 2001/02/15 02:12:22 julian Exp $";
#endif lint
#include <string.h>
#include "compat.h"
#include "f77types.h"

	double	atoangle PARMS((char*));

/*
 * f77-callable interface to atoangle function
 */
float
chtof_(p)
const	char	*p;	/* Character string		(input)	*/
{
	return atof(p);
}
