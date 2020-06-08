#include "SpeciesIO.h"

SpeciesIO::SpeciesIO(){
// Initializes class
//
// input: none
// output:none

	message = "Class var";

}

void SpeciesIO::setMessage(string value){
// Sets attribute hello world
//
// input: string value
// output:none

	message = value;

}

string SpeciesIO::getMessage(){
// Returns hello world
//
// input: none
// output: string value

	return message;
}
