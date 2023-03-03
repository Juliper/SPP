
/*
Double homogeneous model ray tracer, Copyright (c) 2001-2002, University of Amsterdam

Copying, use and development for education and research purposes 
permitted as long as this license is not removed from the files.

All rights for commercial use reserved; for more information
contact Daniel Fontijne (fontijne@science.uva.nl)

This software is unsupported.
*/

#ifndef _U_TIME_H_
#define _U_TIME_H_

#ifdef __cplusplus
extern "C" {
#endif

#define UTERR_OK -1

/*
returns time (in seconds) passed, since some arbitrary point in the past (e.g. client or workstation startup)
returns -1.0 on error
*/
double u_timeGet();

/*
returns time (in seconds) passed since last call to this function (returns 0.0 the first time)
returns -1.0 on error
*/
double u_timeElapsed();

#ifdef __cplusplus
}
#endif

#endif /* _U_TIME_H_ */

