SexesProceedings
================

This programme is written in C. Each pair of .c and .h file is written for different parts of the life cycle. The .c files include the relevant code and the .h files the relevant headers. All functions are called in the main.c (and MainLarge.c for the model with mating type mutants) and you may run the program by compiling all .c files. Below is a short description of each file.

*
*
*

BiogenesisFunctions.c and BiogenesisFunctiosn.h
This includes probability functions that can be used to model selfish mutations. We assume that mitochondria undergo some form of biogenesis where mutants have a replicative advantage. The Biogenesis function is called in the main. The .h file includes all function headers.

EnvironmentalFluctuations.c and EnvironmentalFluctuations.h
This includes functions that are used to force the population to switch the optimal mitonuclear state periodically. Fitness depends on whether the mitonuclear state matches an environmental element and this changes periodically. There are two versions of most functions, one to be used in the simple case without mating types and one to be used with mating types (all functions that include Large in their name). The .h file includes all function headers.
 
GeneralFunctions.c and GeneralFunctions.h
A collection of simple mathematical functions used throughout the program. The .h file includes all function headers.

Global.h
This is a header file where the model parameters are defined. 

InitiationFunctions.c and InitiationFunctions.h
This includes initiations of all matrices used. The .h file includes all function headers.

Main.c
This is the main file for the simple model without mating types. All functions are called in the right sequence according to the life cycle. 

MainLarge.c
This is the main file for the more complex model with mating types. All functions are called in the right sequence according to the life cycle. Different scenarios of mutations etc. can be considered by including selection of the genotype matrices in the programme.

MeiosisFunctions.c and MeiosisFunctions.h
This includes functions that incorporate the probability transformations to implement meiosis. The .h file includes all function headers.

MutationsFunctions.c and MutationsFunctions.h
This includes functions that incorporate the probability transformations to implement mutation.  The .h file includes all function headers.

SelectionFunctions.c and SelectionFunctions.h
This includes functions that incorporate the probability transformations to implement selection. All the functions that include Large in their names are used for the more complex model with mating types (MainLarge.c).  The .h file includes all function headers.

SyngamyFunctions.c and SyngamyFunctions.h    
This includes functions that incorporate the probability matrices for syngamy. Syngamy1 is for uniparental inheritance and syngmay2 for biparental inheritance. The .h file includes all function headers.


