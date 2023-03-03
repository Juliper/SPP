#include <string.h>
#include <string>

std::string to_string(int i) {
  char buf[32];
  sprintf(buf, "%d", i);
  return buf;
}

std::string to_string(double d) {
  char buf[64];
  sprintf(buf, "%2.20e", d);
  return buf;
}

std::string to_string(float f) {
  char buf[64];
  sprintf(buf, "%2.20e", f);
  return buf;
}

/*
reads line up to and including newline character
returns line _excluding_ newline character in 'buf'
return length of line, or negative on failure
removes crap from end of line
returns number of tabs a start of line
removes crap from start of line
'n' is the maximum length of this line
*/
int read_line(FILE* F, char* buf, int n, int* ntabs) {
  int i;
  if (fgets(buf, n - 1, F) == 0) {
    buf[0] = 0;
    return -1;
  }
  n = (int)strlen(buf);
  while ((n > 0) && (buf[n - 1] <= 32))
    n--;
  buf[n] = 0;

  *ntabs = i = 0;
  while ((buf[i] <= 32) && (i < n)) {
    if (buf[i] == '\t')
      (*ntabs)++;
    i++;
  }

  memmove(buf, buf + i, n - i + 1);

  return (n) ? (n - 1) : 0;
}

/*
Returns true if 'line' is an empty line
*/
bool empty_line(const char* line) {
  return ((line[0] == '#') || (line[0] == '\0'));
}

int read_filename(char* line, char* s1) {
  char entryname[256];
  char *p1, *p2;
  if (p1 = strchr(line, '\"')) {
    if (p2 = strchr(p1 + 1, '\"')) {
      memcpy(s1, p1 + 1, p2 - p1 - 1);
      s1[p2 - p1 - 1] = 0;
    } else {
      return -1;
    }
  } else {
    if (sscanf(line, "%s %s", entryname, s1) != 2) {
      return -1;
    }
  }
  return 0;
}

std::string getMVString(const std::string& str, char startDL, char endDL) {
  int startIdx = (int)str.find_first_of(startDL);
  int endIdx   = (int)str.find_last_of(endDL);
  if ((startIdx == str.npos) || (endIdx == str.npos))
    return "";
  else
    return str.substr(startIdx + 1, endIdx - startIdx - 1);
}
