#ifndef UTIL_H
#define UTIL_H

#include <stdarg.h>
#include <stdio.h>

// modify this if needed
#define UTIL_H_DEBUG
#define NRNDET_DEBUG


#ifdef UTIL_H_DEBUG
	#define PRINT_DEBUG(content) \
		fprintf(stderr, content);
	#define PRINTF_DEBUG(content, ...) \
		fprintf(stderr, content, __VA_ARGS__);
#else
	#define PRINT_DEBUG(content)
	#define PRINTF_DEBUG(content, ...)
#endif


#ifdef NRNDET_DEBUG
	#define PRINTF_DEBUG_NRNDET(content, ...)
#endif

#endif