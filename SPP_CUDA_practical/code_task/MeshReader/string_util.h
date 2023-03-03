#ifndef STRING_UTIL_H_
#define STRING_UTIL_H_
//#define WIN32

#include <string>
#include <string.h>

std::string to_string(int i);
std::string to_string(double d);
std::string to_string(float f);

// string compare without case
inline int strcmpnc(const char *s1, const char *s2) {
#ifdef WIN32
	return _stricmp(s1, s2);
#else // the rest of the world:
	return strcasecmp(s1, s2);
#endif
}

int read_line(FILE *F, char *buf, int n, int *ntabs);
bool empty_line(const char *line);
int read_filename(char *line, char *s1);

std::string getMVString(const std::string &str, char startDL, char endDL);

#endif /* STRING_UTIL_H_ */

