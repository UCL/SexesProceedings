#ifndef GLOBAL
#define GLOBAL //Global.h

// Dimensions of the tables
#define NUC (3)                 //  number of nuclear states	
#define MIT (50)				// number of mitochondria

//my parameters
#define mu (0.01)				//miotochondrial mutation rate 0 to 1
#define mub (0.0)				//miotochondrial mutation rate 1 to 0
#define nu (0.000)				//nuclear mutation rate 0 to 1
#define nub (0.000)				//nuclear mutation rate 1 to 0
#define runs (2000)				//number of generations 
#define M (50)
#define types (3)				//number of genotypes
#define E (0.75)					//cost to be implemented duting environmental fluctuations
#define changenuc (200)			//every how many generations to switch states. 
#define duration (21)			//duration of cost implementation
#define DivisionNum (1)			//number of division steps that development needs
#define genotypes (8)			//number of genotype in run (3 for small, 8 for large)


//Output file
#define MEMORYFILE "memory.txt" //to write ratio results
#define FITNESSFILE "mean.txt" //to write averages
#define VARIANCEFILE "variance.txt" //to write averages
#define aafile "a.txt" //to write results
#define AAfile "b.txt" //to write results
#define Aafile "c.txt" //to write results
#define homoz file "hom.txt" //to write results

#endif
