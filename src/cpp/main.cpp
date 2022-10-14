//============================================================================
// Name        : main.cpp
// Author      : Mario Levorato
// Version     :
// Copyright   : Copyright (c) 2013
// Description : Program entry point. Includes main function and CLI.
//============================================================================

#include "controller/include/CommandLineInterfaceController.h"

using namespace controller;

int main(int ac, char* av[])
{
    CommandLineInterfaceController controller;
	int return_value = controller.processArgumentsAndExecute(ac, av);

	return return_value;
}
