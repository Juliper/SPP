#include "Log.h"

#include <stdarg.h>

Log* Log::_instance = 0;

Log* Log::Instance() {
  if (_instance == 0) {
    _instance = new Log("test.txt");
  }
  return _instance;
}

Log* Log::Instance(std::string filename) {
  if (_instance == 0) {
    _instance = new Log(filename.c_str());
  }
  return _instance;
}

Log::Log(std::string filename) {
  m_stream.open(filename.c_str(), ios::out | ios::app);
}

void Log::Write(const char* logline, ...) {
  va_list argList;
  char cbuffer[1024];
  va_start(argList, logline);
  vsnprintf(cbuffer, 1024, logline, argList);
  va_end(argList);
  m_stream << cbuffer << endl;
}

void Log::Clear() {
  m_stream.close();
  m_stream.open("test.txt", ios::out | ios::trunc);
}

Log::~Log() {
  m_stream.close();
}