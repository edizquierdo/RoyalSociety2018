// Integration parameters
const double Duration = 20.0;           // Seconds
const double Transient = 10.0;          //
const double StepSize = 0.01;
const int N_curvs = 23;                 // Number of cuvature points

// Used for Dumping: Frame rate for recording datais set to 50 frames per second
const double fps = 25.0;
const int skip = (int) (1/(StepSize*fps));

// Genotype -> Phenotype Mapping (Ventral cord)
const double	  BiasRange			         	= 15.0;
const double    SCRange                 = 15.0;
const double    CSRange                 = 15.0;
const double    TauMin                  = 0.5; //
const double    TauMax                  = 2.0;

const double    ESRange                 = 2.0;

const double    SRmax                   = 200.0;
const double    NMJmax                  = 1.0;

// (Head)
const double    HCSRange                = 15.0;

// Fitness
const double    AvgSpeed = 0.00022;              // Average speed of the worm in meters per seconds
const double    BBCfit = AvgSpeed*Duration;

// Size of genotype (VC)
int	VectSize = 30;