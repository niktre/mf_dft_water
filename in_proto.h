void AllocArrays (void);

void AllocSubstrate (void);

void CalcFreeEn (int, int, int);

void CalcPartSumQ (int);

void CheckConverge (int);

void CheckSuccess (double, double, double);

void DefineSimpson (void);

void DivideSpace (void);

void InitParameters (void);

void PassConsoleParams (int, char **);

void PrintSnapField (int, int);

void reverse(char *);

void StoreParameters (void);

void CommHaloYZ (void);
void ReadSubFile (void);
void ReadSubFunc (void);
void PrintSubFile (void);

/* STRESS TENSOR FUNCTIONS */
void CalcForceZ (void);
void CalcSigmas (void);
void CalcDivSigma (void);
void PrintDiff (void);