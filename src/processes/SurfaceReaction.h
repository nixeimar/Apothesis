#include <string>
using namespace std;

class SurfaceReaction
{
	public:
		SurfaceReaction();
		string getMessage();
		void setMessage(string);
		
	private:
		string message;
};

