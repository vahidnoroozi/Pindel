#ifndef __MESSAGE_H
#define __MESSAGE_H

#include <stdio.h>
#include <stdarg.h>

enum {NO_MSG, INFO_MSG, DEBUG_MSG, WARNING_MSG, ERROR_MSG};
enum {NO_ERROR, CUSTOM_ERROR, NO_DATA, MEMORY_ALLOCATION, FILE_NOT_FOUND,
	FILE_OP_FAILED, END_OF_FILE, FILE_FORMAT_ERROR, INVALID_CMD_OPTION,
	INVALID_CMD_ARGUMENT, RESOURCE_NOT_FOUND, INTERNAL_MISMATCH, 
	INTERNAL_ERROR, STATE_SPACE_OVERFLOW};
enum {SILENT, QUIET, MINIMAL, VERBOSE};

#ifdef INFO_WITH_LINO
#define INFO(...) message(stderr, __FILE__, __LINE__, INFO_MSG, \
		NO_ERROR, __VA_ARGS__)
#else
#define INFO(...) message(stderr, NULL, 0, INFO_MSG, \
		NO_ERROR, __VA_ARGS__)
#endif

int message(FILE *fp, const char *file_name, int line,
	int msg_type, int msg_id, const char *format, ...);

#endif
