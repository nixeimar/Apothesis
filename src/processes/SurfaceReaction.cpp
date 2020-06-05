#include "SurfaceReaction.h"

SurfaceReaction::SurfaceReaction(){
// Initializes class
//
// input: none
// output:none

	message = "Class var";

}

void SurfaceReaction::setMessage(string value){
// Sets attribute hello world
//
// input: string value
// output:none

	message = value;

}

string SurfaceReaction::getMessage(){
// Returns hello world
//
// input: none
// output: string value

	return message;
}
