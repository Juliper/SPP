#ifndef LOG_H
#define LOG_H

#include <fstream>
#include <string>
#include <cstring>

using namespace std;

class Log {
public:
	static Log* Instance();
	static Log* Instance(string filename); 
	Log(std::string filename);
	~Log();
	void Write(const char* logline,...);
	void Clear();

private:
	static Log* _instance;
	ofstream m_stream;
};

#endif