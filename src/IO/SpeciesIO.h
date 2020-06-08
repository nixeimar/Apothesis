#include <string>
#include "io.h"

using namespace std;

class IO;

class SpeciesIO
{
	public:
		SpeciesIO();
		string getMessage();
		void setMessage(string);
		
	private:
		string message;
};

