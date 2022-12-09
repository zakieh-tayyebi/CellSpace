#include<iostream>
#include<string>
using namespace std;

int main(int argc, char** argv){
	string line, barcode = string(argv[1]);
	while(getline(cin, line)) cout << line << "\t" << "CB:Z:" << barcode << endl;
	
	return 0;
}
