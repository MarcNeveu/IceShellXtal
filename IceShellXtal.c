/*
 * IceShellXtal.c
 *
 *  Created on: May 31, 2022
 *      Author: mneveu
 */

#include <unistd.h>    // To check current working directory
#include <math.h>
#include <IPhreeqc.h>  // To use the external PHREEQC geochemical code
#include <Var.h>       // To use the external PHREEQC geochemical code

#include "Plot.h"

int UpdateDisplays (SDL_Renderer* renderer, SDL_Texture* background_tex, SDL_Texture* pies_tex, SDL_Texture **num_tex, char* FontFile);
int ChamberPlot(SDL_Surface **chamber, char *FontFile, int ntemp, double radius, double **R2Hcompo, double Rchamber, double R1, int nsalts, int contour, int dbase_type);
int ExtractWrite(int instance, double*** data, int line, int nvar);
int MineralName(int i, char **name, int dbase_type);
double Vm(char *name, char *dbase);
double VsphSegm(double x, double R2, double b);
double R2iceCap (double x, double Vsol, double R1, double H);
SDL_Color minColor (int i, int dbase_type);

int main(int argc, char *argv[]){

	int i = 0;
	int j = 0;

	//-------------------------------------------------------------------
	// Startup
	//-------------------------------------------------------------------

	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("IceShellXtal v22.6\n");
	printf("-------------------------------------------------------------------\n");

	// Get current directory. Works for Mac only! To switch between platforms, see:
	// http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe

	char path[1024];
	unsigned int size = sizeof(path);
	path[0] = '\0';

	if (_NSGetExecutablePath(path, &size) == 0)
		printf("\n");
	else
	    printf("IceShellXtal: Couldn't retrieve executable directory. Buffer too small; need size %u\n", size);

	//-------------------------------------------------------------------
	// Read inputs
	//-------------------------------------------------------------------
	printf("\n------------------------------------------------------\n");
	printf("Reading input file...\n");
	printf("------------------------------------------------------\n\n");

	char infile[1024];    infile[0] = '\0';    // Path to input file
	int dbase_type = 0;     // 0 = TC17 (chlorides and sulfates only), 1 = also Carbonates, Si, CH4 TODO make this an input
	int ntemp = 0; 		    // Number of temperature steps in simulation TODO make this an input, set to last successful i + 1
	double temp = 0.0;      // Simulation temperature (Celsius)
	double temp_step = 0.0; // Temperature step (C or K)

	int eqFrac = 0;         // 0 = equilibrium crystallization, 1 = fractional crystallization
	int contour = 0;        // Plot salt deposit contours?

	FILE *f;
	char *PHRQin = (char *) malloc(1024*sizeof(char));
	PHRQin[0] = '\0';
	strncat(infile,path,strlen(path)-20);
	strcat(infile,"IceShellXtal_input.txt");

	f = fopen (infile,"r"); 	// Open input file
	if (f == NULL) printf("IceShellXtal: Missing input file. Path: %s\n", infile);
	int tab = 44;
	fseek(f,196,SEEK_SET);
	fseek(f,tab,SEEK_CUR); fscanf(f, "%s", PHRQin);      printf("PHREEQC input file                        | %s\n", PHRQin);
	fseek(f,tab,SEEK_CUR); fscanf(f, "%d", &dbase_type); printf("PHREEQC database: 0=TC17, 1=frezchemSiCH4 | %d\n", dbase_type);
	fseek(f,tab,SEEK_CUR); fscanf(f, "%d", &ntemp);      printf("Number of temperature steps               | %d\n", ntemp);
	fseek(f,tab,SEEK_CUR); fscanf(f, "%lg", &temp);      printf("Starting temperature (C)                  | %g\n", temp);
	fseek(f,tab,SEEK_CUR); fscanf(f, "%lg", &temp_step); printf("Temperature step, < 0                     | %g\n", temp_step);
	fseek(f,tab,SEEK_CUR); fscanf(f, "%d", &eqFrac);     printf("Crystallization mode: 0=equil, 1=frac     | %d\n", eqFrac);
	fseek(f,tab,SEEK_CUR); fscanf(f, "%d", &contour);    printf("Plot salt deposit contours? 0=no, 1=yes   | %d\n", contour);

	fclose(f);

	infile[0] = '\0';
	strncat(infile,path,strlen(path)-20);
	strcat(infile, PHRQin);

	free (PHRQin);

	//-------------------------------------------------------------------
	// Run PHREEQC
	//-------------------------------------------------------------------
	printf("\n------------------------------------------------------\n");
	printf("Running PHREEQC...\n");
	printf("------------------------------------------------------\n\n");

	int phreeqc = 0;
	char dbase[1024];     dbase[0] = '\0'; 	// Path to thermodynamic database
	int nvar = 200; // Max. number of variables output by PHREEQC SELECTED_OUTPUT

	double **simdata = (double**) malloc(ntemp*sizeof(double*)); // Simulation data storage table
	if (simdata == NULL) printf("IceShellXtal: Not enough memory to create simdata[ntemp]\n");
	for (i=0;i<ntemp;i++) {
		simdata[i] = (double*) malloc(nvar*sizeof(double));
		if (simdata[i] == NULL) printf("IceShellXtal: Not enough memory to create simdata[ntemp][nvar]\n");
	}
	for (i=0;i<ntemp;i++) {
		for (j=0;j<nvar;j++) simdata[i][j] = 0.0;
	}

	char temperature[20];

	strncat(dbase,path,strlen(path)-20);
	if (!dbase_type) strcat(dbase,"PHREEQC/PHREEQC_TC17_vm.txt");
	else strcat(dbase,"PHREEQC/frezchemSiCH4.dat");

	phreeqc = CreateIPhreeqc(); // Run PHREEQC
	if (LoadDatabase(phreeqc,dbase) != 0) OutputErrorString(phreeqc);
	SetSelectedOutputFileOn(phreeqc,1);
	SetDumpStringOn(phreeqc, 1);
	SetDumpFileOn(phreeqc, 1);
	printf("Running PHREEQC\n"); // First temperature step is run from the PHREEQC input file
	printf("%s\n", infile);
	if (RunFile(phreeqc,infile) != 0) OutputErrorString(phreeqc);
	else printf("PHREEQC ran successfully\n");

	ExtractWrite(phreeqc, &simdata, 1, nvar); // Extract SELECTED_OUTPUT into simdata

	infile[0] = '\0';
	strncat(infile,path,strlen(path)-20);
	strcat(infile, "Output_compo.txt");
	f = fopen (infile,"w"); // Open output file
	// Print headers to output file
	fprintf(f, "Sim \tpH \tT(C) \tIS \tmH2O(kg) \tchrgbal \tpcterr \tCa(M) \tMg(M) \tNa(M) \tK(M) \tCl(M) \tS(M) \t");
	if (dbase_type) { // Carbonate-Si-CH4 database
		fprintf(f, "C(M) \tSi(M) \tMtg(M) \t");
	}
	fprintf(f, "Bah(M) \tlog_aH2O \t"
				"Ice(mol) \tdIce(mol) \tHalite(mol) \tdHalite(mol) \tHydrohalite(mol) \tdHydrohalite(mol) \tSylvite(mol) \tdSylvite(mol) "
				"\tAntarcticite(mol) \tdAntarcticite(mol) \tBischofite(mol) \tdBischofite(mol) \tMgCl2:8H2O(mol) \tdMgCl2:8H2O(mol) \t"
				"MgCl2:12H2O(mol) \tdMgCl2:12H2O(mol) \tCarnallite(mol) \tdCarnallite(mol) \tTachyhydrite(mol) \tdTachyhydrite(mol) \t"
				"Anhydrite(mol) \tdAnhydrite(mol) \tArcanite(mol) \tdArcanite(mol) \tEpsomite(mol) \tdEpsomite(mol) \tGypsum(mol) \tdGypsum(mol) \t"
				"Meridianiite(mol) \tdMeridianiite(mol) \tMirabilite(mol) \tdMirabilite(mol) \tPicromerite(mol) \tdPicromerite(mol) \t"
				"Glaserite(mol) \tdGlaserite(mol) \tBloedite(mol) \tdBloedite(mol) \tKieserite(mol) \tdKieserite(mol) \tThenardite(mol) \t"
				"dThenardite(mol) \tHexahydrite(mol) \tdHexahydrite(mol) \t");
	if (!dbase_type) { // TC17 database
		fprintf(f, "Na2SO4:7H2O(mol) \tdNa2SO4:7H2O(mol) \tStarkeyite(mol) \tdStarkeyite(mol) \tPentahydrite(mol) \tdPentahydrite(mol) \t"
				"Bassanite(mol) \tdBassanite(mol) \tGlauberite(mol) \tdGlauberite(mol) \tLabile_Salt(mol) \tdLabile_Salt(mol) \tLeonite(mol) \tdLeonite(mol) \t"
				"Langbeinite(mol) \tdLangbeinite(mol) \tSyngenite(mol) \tdSyngenite(mol) \tPolyhalite(mol) \tdPolyhalite(mol) \tKainite(mol) \tdKainite(mol) \t");
	}
	else { // Carbonate-Si-CH4 database
		fprintf(f, "Aragonite(mol) \tdAragonite(mol) \tArtinite(mol) \tdArtinite(mol) \tCalcite(mol) \tdCalcite(mol) \tDolomite(mol) \tdDolomite(mol) \t"
				"Hydromagnesite(mol) \tdHydromagnesite(mol) \tIkaite(mol) \tdIkaite(mol) \tKalicinite(mol) \tdKalicinite(mol) \tLandsfordite(mol) \t"
				"dLandsfordite(mol) \tMagnesite(mol) \tdMagnesite(mol) \tNa2CO3:7H2O(mol) \tdNa2CO3:7H2O(mol) \tNahcolite(mol) \tdNahcolite(mol) \t"
				"Natron(mol) \tdNatron(mol) \tNesquehonite(mol) \tdNesquehonite(mol) \tTrona(mol) \tdTrona(mol) \tVaterite(mol) \tdVaterite(mol) \t"
				"Akermanite(mol) \tdAkermanite(mol) \tAnthophyllite(mol) \tdAnthophyllite(mol) \tAntigorite(mol) \tdAntigorite(mol) \tChalcedony(mol) \t"
				"dChalcedony(mol) \tChrysotile(mol) \tdChrysotile(mol) \tDiopside(mol) \tdDiopside(mol) \tEnstatite(mol) \tdEnstatite(mol) \tForsterite(mol) \t"
				"dForsterite(mol) \tQuartz(mol) \tdQuartz(mol) \tSepiolite(mol) \tdSepiolite(mol) \tSiO2(am)(mol) \tdSiO2(am)(mol) \tTalc(mol) \t"
				"dTalc(mol) \tMethane_hydrate(mol) \tdMethane_hydrate(mol) \t");
	}
	fprintf(f, "SolnDensity(kg/L) \tSolnVolume(L) \tPressure(atm) \tConductivity(uS/cm)\n");
	for (i=0;i<nvar;i++) {
		fprintf(f, "%g\t", simdata[1][i]); // Print to output file
		printf("%g\t", simdata[1][i]);     // Print to terminal
	}
	fprintf(f, "\n");
	printf("\n");

	for (i=2;i<ntemp;i++) { // Re-run the same input at lower temperatures, using the solution of the previous temperature step.
		AccumulateLine(phreeqc, "TITLE Next temperature iteration\n");
		AccumulateLine(phreeqc, "SOLUTION_MODIFY 1");
		temp += temp_step;
		sprintf(temperature, "\ttemp %g", temp);
		AccumulateLine(phreeqc, temperature);
		AccumulateLine(phreeqc, "USE solution 1");
		AccumulateLine(phreeqc, "USE equilibrium_phases 1"); // Allowing the same phases to precipitate. Amounts of precipitate are memorized only in equilibrium crystallization mode.

		AccumulateLine(phreeqc, "SELECTED_OUTPUT"); // The SELECTED_OUTPUT instructions below should be exactly the same as those in the PHREEQC input file, otherwise the first row of outputs will be different.
		AccumulateLine(phreeqc, " -file ./selected.out");
		AccumulateLine(phreeqc, " -high_precision true");
		AccumulateLine(phreeqc, " -reset false");
		AccumulateLine(phreeqc, " -simulation true");
		AccumulateLine(phreeqc, " -pH true");
		AccumulateLine(phreeqc, " -temperature true");
		AccumulateLine(phreeqc, " -pe false");
		AccumulateLine(phreeqc, " -ionic_strength true");
		AccumulateLine(phreeqc, " -water true");
		AccumulateLine(phreeqc, " -charge_balance true");
		AccumulateLine(phreeqc, " -percent_error true");
		if (!dbase_type) AccumulateLine(phreeqc, " -totals Ca Mg Na K Cl S Bah");
		else AccumulateLine(phreeqc, " -totals Ca Mg Na K Cl S C Si Mtg Bah");
		AccumulateLine(phreeqc, " -activities H2O");
		if (!dbase_type) // TC17 database
			AccumulateLine(phreeqc, " -equilibrium_phases Ice(s) Halite Hydrohalite Sylvite Antarcticite Bischofite MgCl2:8H2O "
				"MgCl2:12H2O Carnallite Tachyhydrite Anhydrite Arcanite Epsomite Gypsum Meridianiite Mirabilite Picromerite "
				"Glaserite Bloedite Kieserite Thenardite Hexahydrite Na2SO4:7H2O Starkeyite Pentahydrite Bassanite "
				"Glauberite Labile_Salt Leonite Langbeinite Syngenite Polyhalite Kainite");
		else // carbonate-Si-CH4 database
			AccumulateLine(phreeqc, " -equilibrium_phases Ice(s) Halite Hydrohalite Sylvite Antarcticite Bischofite MgCl2:8H2O "
							"MgCl2:12H2O Carnallite Tachyhydrite Anhydrite Arcanite Epsomite Gypsum Meridianiite Mirabilite Picromerite "
							"Glaserite Bloedite Kieserite Thenardite Hexahydrite Aragonite Artinite Calcite Dolomite Hydromagnesite "
							"Ikaite Kalicinite Landsfordite Magnesite Na2CO3:7H2O Nahcolite Natron Nesquehonite Trona Vaterite Akermanite "
							"Anthophyllite Antigorite Chalcedony Chrysotile Diopside Enstatite Forsterite Quartz Sepiolite SiO2(am) Talc Methane_hydrate");

		AccumulateLine(phreeqc, "USER_PUNCH");
		AccumulateLine(phreeqc, " -headings Solution_density Solution_volume Pressure Specific_conductance");
		AccumulateLine(phreeqc, " -start");
		AccumulateLine(phreeqc, " PUNCH RHO");      // kg/L
		AccumulateLine(phreeqc, " PUNCH SOLN_VOL"); // L
		AccumulateLine(phreeqc, " PUNCH PRESSURE"); // atm
		AccumulateLine(phreeqc, " PUNCH SC");       // uS/cm, need to have -dw values in database for nonzero output

		AccumulateLine(phreeqc, "SAVE solution 1");
		if (!eqFrac) AccumulateLine(phreeqc, "SAVE equilibrium_phases 1"); // Save solids for equilibrium fractionation
		AccumulateLine(phreeqc, "END");

		printf("Running PHREEQC at T=%g C\n", temp);
		if (RunAccumulated(phreeqc) != 0) OutputErrorString(phreeqc);
		else printf("PHREEQC ran successfully at i=%d (ntemp=%d)\n", i, ntemp);

		ExtractWrite(phreeqc, &simdata, i, nvar); // Memorize PHREEQC SELECTED_OUTPUT in simdata
		for (j=0;j<nvar;j++) {
			fprintf(f, "%g\t", simdata[i][j]); // Print to output file
			printf("%g\t", simdata[i][j]);     // Print to terminal
		}
		fprintf(f, "\n");
		printf("\n");
	}
	fclose(f);

	if (DestroyIPhreeqc(phreeqc) != IPQ_OK) OutputErrorString(phreeqc);

	//-------------------------------------------------------------------
	// Analyze PHREEQC output
	//-------------------------------------------------------------------
	printf("\n------------------------------------------------------\n");
	printf("Analyzing PHREEQC output...\n");
	printf("------------------------------------------------------\n\n");

	// Find salts that formed
	int mineral_index[nvar];                       // Index of mineral
	int k = 0;                                     // Salt counter
	char *name = (char *)malloc(200*sizeof(char)); // Name of mineral

	for (k=0;k<nvar;k++) mineral_index[k] = 0;
	k = 0;
	name[0] = '\0';

	int isaltmin = 0; // Min column index of salts in PHREEQC output
	int isaltmax = 0; // Max column index of salts in PHREEQC output

	if (!dbase_type) { // ... in case of TC17 database
		isaltmin = 15;
		isaltmax = isaltmin+2*33;
	}
	else {             // ... for carbonate-Si-CH4 database
		isaltmin = 18;
		isaltmax = isaltmin+2*50;
	}

	// Find number and indices of minerals that crystallized
	for (j=isaltmin;j<isaltmax;j=j+2) { // Only for solids columns, skip delta columns
		for (i=1;i<ntemp;i++) simdata[0][j] += simdata[i][j]; // Sum in first row
		if (simdata[0][j] > 0.0) {
			mineral_index[k] = j;
			k++;
		}
	}
	int nsalts = k-1;
	printf("%d solids crystallized:\n", nsalts+1);
	for (i=0;i<nsalts+1;i++) {
		MineralName(mineral_index[i], &name, dbase_type);
		printf("%s, index %d, %g mol, Vm=%g cm3/mol\n", name, mineral_index[i], simdata[0][mineral_index[i]], Vm(name, dbase));
	}

	// Convert moles of salts formed to volumes
	// Build table of temperatures, total volume, solution volume, total salt volume, ice volume, and individual salt volumes at each temperature
	double **volumes = (double**) malloc(ntemp*sizeof(double*));
	if (volumes == NULL) printf("IceShellXtal: Not enough memory to create volumes[%d]\n", ntemp);
	for (i=0;i<ntemp;i++) {
		volumes[i] = (double*) malloc((nsalts+5)*sizeof(double));
		if (volumes[i] == NULL) printf("IceShellXtal: Not enough memory to create volumes[%d][%d]\n", ntemp, nsalts+5); // 18 columns in System_main_tbl.txt
	}
	for (i=0;i<ntemp;i++) {
		for(j=0;j<nsalts+5;j++) volumes[i][j] = 0.0;
	}
	for (i=1;i<ntemp;i++) {
		volumes[i][0] = simdata[i][2] + 273.15; // Temperature (K)
		volumes[i][2] = simdata[i][isaltmax+1] *1000.0; // Solution volume (cm3)
		for (j=0;j<nsalts+1;j++) {
			MineralName(mineral_index[j], &name, dbase_type);
			volumes[i][j+4] = simdata[i][mineral_index[j]]*Vm(name, dbase); // Individual solid volumes (cm3)
			if (j>0) volumes[i][3] += volumes[i][j+4];   // Total salt volume (cm3), excluding j=0 which is Ice(s) per SELECTED_OUTPUT instructions
		}
		volumes[i][1] = volumes[i-1][1] - volumes[i-1][2] + volumes[i][2] + volumes[i][3] + volumes[i][4]; // Total chamber volume (cm3) = old vtot - old_solution + solution + total new salt + new ice
	}

	if (!eqFrac) { // Equilibrium well-mixed mode
		printf("\nT(K) \tVsol \tV_ice \tIndividual salt volumes (cm3)\n");
		for (i=0;i<ntemp;i++) {
			printf("%g \t", volumes[i][0]);
			printf("%g \t", volumes[i][2]);
			for(j=0;j<nsalts+1;j++) {
				printf("%g \t", volumes[i][j+4]);
			}
			printf("\n");
		}
	}
	else { // Fractional, no mixing mode
		printf("\nT(K) \tVtot \tVsol \tVsalt \tV_ice \tIndividual salt volumes (cm3)\n");
		for (i=0;i<ntemp;i++) {
			for(j=0;j<nsalts+5;j++) {
				printf("%g \t", volumes[i][j]);
			}
			printf("\n");
		}
	}

	// Build table of temperatures, R2/R1, H/R1, and vol% salt compositions at each temperature
	double **R2Hcompo = (double**) malloc(ntemp*sizeof(double*));
	if (R2Hcompo == NULL) printf("IceShellXtal: Not enough memory to create R2Hcompo[%d]\n", ntemp);
	for (i=0;i<ntemp;i++) {
		R2Hcompo[i] = (double*) malloc((nsalts+3)*sizeof(double));
		if (R2Hcompo[i] == NULL) printf("IceShellXtal: Not enough memory to create R2Hcompo[%d][%d]\n", ntemp, nsalts+3); // 18 columns in System_main_tbl.txt
	}
	for (i=0;i<ntemp;i++) {
		for(j=0;j<nsalts+3;j++) R2Hcompo[i][j] = 0.0;
	}
	R2Hcompo[0][1] = 1.0; // = R2/R1(init) and R2Hcompo[0][2] = 0.0 = H/R1(init)
	for (j=0;j<nsalts;j++) R2Hcompo[0][j+3] = mineral_index[j+1];

	for (i=1;i<ntemp;i++) {
		R2Hcompo[i][0] = volumes[i][0]; // Temperature (K)
		for (j=0;j<nsalts;j++) {
			if (volumes[i][3] > 0.0) R2Hcompo[i][j+3] = volumes[i][j+5]/volumes[i][3]*100.0; // %vol of individual salts relative to total salt volume, excluding ice
			else R2Hcompo[i][j+3] = 0.0;
		}
	}

	// Calculate chamber dimensions
	double Rchamber = 0.0; // Chamber outer radius (cm3)
	double R1 = 0.0;       // Inner radius (cm) of ice ring at the temperature of first salt formation (cm3) (fractional not mixed mode) or central pocket at end of computation (equilibrium well-mixed mode)
	int firstSalt = 0;     // Index at which salt first crystallizes

	Rchamber = pow(volumes[1][2]*3.0/4.0/M_PI, 1.0/3.0); // Initial solution volume
	for (i=1;i<ntemp;i++) {
		if (volumes[i][3] > 0.0 && volumes[i][4] > 0.0) { // Both new salt and new ice need to be positive; if new salt happens first, that salt drains to the ocean
			firstSalt = i;
			R1 = pow(volumes[i-1][2]*3.0/4.0/M_PI, 1.0/3.0); // Solution volume just before the first salts form. Neglect variation in R1 with freezing; otw, will complicate resets of R1 if R2 is ever < R1-H
			break;
		}
	}

	if (!eqFrac) { // Well-mixed, equilibrium crystallization mode. Output central pocket composition and stop here.
		double VctrPocket = 0.0; // Central pocket volume (cm3)
		for (j=0;j<nsalts;j++) VctrPocket += volumes[ntemp-1][j+5]; // Total cumulative salt volume (cm3)
		VctrPocket += volumes[ntemp-1][2];                          // Solution volume at last temperature step (cm3)
		R1 = pow(VctrPocket*3.0/4.0/M_PI, 1.0/3.0);                 // Radius of central pocket (cm)

		printf("\nWell-mixed, equilibrium crystallization mode.\n"
				"Chamber of radius Rchamber is pure ice outside R1.\n"
				"Composition of central pocket inside R1 (vol%%):\n");
		printf("Residual briny solution: %g\n", volumes[ntemp-1][2]/VctrPocket*100.0);
		for (j=0;j<nsalts;j++) {
			MineralName(mineral_index[j+1], &name, dbase_type);
			printf("%s: %g\n", name, volumes[ntemp-1][j+5]/VctrPocket*100.0);
		}
		printf("Rchamber = %g, R1 = %g\n", Rchamber, R1);

		printf("\nExiting IceShellXtal...\n", Rchamber, R1);
		exit(0);
	}
	// Not mixed, fractional crystallization mode. Output chamber size and inner radius R1 at temperature of first salt formation, and go on.
	printf("\nRchamber = %g, R1 = %g\n", Rchamber, R1);
	printf("Salt first precipitates in the presence of ice at step %d\n", firstSalt);

	// Calculate H and R2
	double deltah = 0.0;          // Incremental salt deposit height (cm)
	double H = 0.0;               // Cumulative salt deposit height (cm)
	double b = 0.0;               // Width of spherical segment (cm) TODO double-check
	double V = 0.0;               // Incremental salt volume (cm3)
	double R2 = R1;               // Radius of remaining solution volume (cm)
	double R1calc = R1;           // Value of R1 used in R2 and H calculation, different from R1 if ice ever covers salt

	double low = 0.0;             // Storage variables for binary search
	double mid = 0.0;
	double high = 0.0;
	double tem = 0.0;
	double Vmid = 0.0;
	double threshold = 1.0e-6;    // Binary search relative threshold to determine convergence
	int iter = 0;                 // Number of iterations in binary search

	for (i=firstSalt;i<ntemp;i++) {
//		printf("i=%d, Vsol=%g, b=%g, H=%g, R1calc=%g, R2=%g\n", i, volumes[i][2], b, H, R1calc, R2);
		// Calculate H, the cumulative height of the salt deposit, using a binary search. Store in R2Hcompo.
	    // If ice covers salt, reset salt deposit height and reset R1calc to R2
	    if (R2 < R1calc-H) {
	        for (j=i;j<ntemp;j++) R1calc = R2;
	        R2 = R1calc;
	        H = 0;
	    }
	    if (volumes[i][2] < volumes[0][2]*1.0e-3) break; // Stop when solution volume is sufficiently low
	    else {
	        b = R1calc - H;
	        V = volumes[i][3]; // Incremental salt volume (cm3)
	        mid = 0.0;
	        tem = 0.0;
	        low = 0.0;
	        high = R2;

	        if (V - VsphSegm(low, R2, b) > 0.0) { // Invert bounds
	        	printf("Inverting bounds in binary search for deltah\n");
	            tem = low;
	            low = high;
	            high = tem;
	        }
	        while (VsphSegm(high, R2, b) < 0.0) { // Decrease high bound
	            high *= (1.0-threshold);
	        }
	        while (VsphSegm(low, R2, b) < 0.0) { // Increase low bound
	            low /= (1.0-threshold);
	        }

	        iter = 0;
	        Vmid = VsphSegm(mid, R2, b);
	        while (fabs(V - Vmid) > threshold*V) {
	            iter++;
	            mid = (low+high)/2.0;
	            Vmid = VsphSegm(mid, R2, b);
	            if(V - Vmid < -threshold*V) low = mid;
	            if (V - Vmid > threshold*V) high = mid;
//	            printf("i=%d, iteration %d, deltah = %.12g cm, V=%g cm3, Vmid = %g cm3\n", i, iter, mid, V, Vmid);
	            if (iter > 100) {
	            	printf("IceShellXtal: could not converge within 100 iterations on H calculation\n");
	            	exit(0);
	            }
	        }
//	        printf("Vmid=%g\n", Vmid);
	    }

	    deltah = mid;          // Incremental salt deposit height (cm)
	    H = H+deltah;          // Cumulative salt deposit height (cm)
	    R2Hcompo[i][2] = H/R1; // Relative to initial R1

	    //Calculate R2, the inner radius of the ice ring that just formed, using a binary search. Store in R2Hcompo.
	    mid = 0.0;
	    tem = 0.0;
	    low = 0.0;
	    high = R2;

		if (R2iceCap(low, volumes[i][2], R1calc, H) > 0.0) {
			tem = low;
			low = high;
			high = tem;
		}

		iter = 0;
		while (fabs(R2iceCap(mid, volumes[i][2], R1calc, H)) > threshold*R2) {
			iter++;
			mid = (low+high)/2.0;
			if(R2iceCap(mid, volumes[i][2], R1calc, H) < -threshold*R2) low = mid;
			if (R2iceCap(mid, volumes[i][2], R1calc, H) > threshold*R2) high = mid;
//			printf("mid=%g R2iceCap=%g\n", mid, R2iceCap(mid, volumes[i][2], R1calc, H));
			if (iter > 100) {
				printf("IceShellXtal: could not converge within 100 iterations on R2 calculation\n");
				exit(0);
			}
		}

	    R2 = mid;
	    R2Hcompo[i][1] = R2/R1; // Relative to initial R1
	}

	// Shift rows of R2Hcompo up to start only when salts first form
	for (i=firstSalt;i<ntemp;i++) {
		for (j=0;j<k+2;j++) R2Hcompo[i-firstSalt+1][j] = R2Hcompo[i][j];
	}
	for (i=ntemp-firstSalt+1;i<ntemp;i++) {
		for (j=0;j<k+2;j++) R2Hcompo[i][j] = 0.0;
	}

	// Print out R2Hcompo
	printf("\nT(K) \tR2/R1 \tH/R1 \tComposition (vol%%)\n");
	for (i=0;i<ntemp;i++) {
		for(j=0;j<nsalts+3;j++) {
			printf("%g \t", R2Hcompo[i][j]);
		}
		printf("\n");
	}

	// Reprint minerals that formed to facilitate identification in printed R2Hcompo table
	printf("%d solids crystallized:\n", nsalts+1);
	for (i=0;i<nsalts+1;i++) {
		MineralName(mineral_index[i], &name, dbase_type);
		printf("%s, index %d, %g mol, Vm=%g cm3/mol\n", name, mineral_index[i], simdata[0][mineral_index[i]], Vm(name, dbase));
	}

	//-------------------------------------------------------------------
	// Launch Sample DirectMedia Layer (SDL) display
	//-------------------------------------------------------------------
	printf("\n------------------------------------------------------\n");
	printf("Plotting chamber...\n");
	printf("------------------------------------------------------\n\n");

	if (SDL_Init(SDL_INIT_EVERYTHING) == -1){
		printf("IceShellXtal: Error: SDL not initialized.");
	}
//	if (TTF_Init() != 0){
//		printf("IceShellXtal: Error: TTF not initialized.");
//	}
	window = SDL_CreateWindow("IceShellXtal", SDL_WINDOWPOS_CENTERED,
		SDL_WINDOWPOS_CENTERED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
	if (window == NULL){
		printf("IceShellXtal: Error: Window not created.");
	}
	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED
		| SDL_RENDERER_PRESENTVSYNC);
	if (renderer == NULL){
		printf("IceShellXtal: Error: Renderer not created.");
	}


	// Display font
	char FontFile[1024];      // Don't forget to free!
	FontFile[0] = '\0';
	strncat(FontFile,path,strlen(path)-20);
	strcat(FontFile,"Graphics/GillSans.ttf");

	//-------------------------------------------------------------------
	//                         Initialize display
	//-------------------------------------------------------------------

    int quit = 0;                       // Close window
	double radius = 200.0;              // Chamber radius (pixels)
	SDL_Texture* background_tex = NULL;
	SDL_Texture* chamber_tex = NULL;
	SDL_Surface* chamber = NULL;

	File2tex("Graphics/IceShellXtalBG/IceShellXtalBG.002.jpeg", &background_tex, path);
	File2surf("Graphics/Transparent.png", &chamber, path);

	SDL_Event e;
	SDL_PollEvent(&e);

	SDL_Texture **num_tex;
	ChamberPlot(&chamber, FontFile, ntemp, radius, R2Hcompo, Rchamber, R1, nsalts, contour, dbase_type);
	chamber_tex = SDL_CreateTextureFromSurface(renderer, chamber);

	//-------------------------------------------------------------------
	//                         Interactive display
	//-------------------------------------------------------------------

	while (!quit){
		while (SDL_PollEvent(&e)){

			if (e.type == SDL_QUIT) quit = 1; // Close window
			if (e.type == SDL_MOUSEBUTTONDOWN) {
				ChamberPlot(&chamber, FontFile, ntemp, radius, R2Hcompo, Rchamber, R1, nsalts, contour, dbase_type);
				chamber_tex = SDL_CreateTextureFromSurface(renderer, chamber);
			}
		}
		UpdateDisplays(renderer, background_tex, chamber_tex, num_tex, FontFile);
	}

	//-------------------------------------------------------------------
	//                      Free remaining mallocs
	//-------------------------------------------------------------------

	SDL_DestroyTexture(background_tex);
	SDL_FreeSurface(chamber);
	SDL_DestroyTexture(chamber_tex);

	//-------------------------------------------------------------------
	// Exit
	//-------------------------------------------------------------------

	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();

	for (i=0;i<ntemp;i++) {
		free(simdata[i]);
		free(volumes[i]);
		free(R2Hcompo[i]);
	}
	free (simdata);
	free(volumes);
	free (R2Hcompo);

	free (name);

	printf("Exiting IceShellXtal...\n");
	return 0;
}

//-------------------------------------------------------------------
//                      Display updating subroutine
//-------------------------------------------------------------------

int UpdateDisplays (SDL_Renderer* renderer, SDL_Texture* background_tex, SDL_Texture* chamber_tex, SDL_Texture **num_tex, char* FontFile) {

	SDL_RenderClear(renderer);
	ApplySurface(0, 0, background_tex, renderer, NULL);
	ApplySurface(0, 0, chamber_tex, renderer, NULL);

//	// Numbers
//	x = 82; y=355; d=52;
//	if (num_tex[j*npressure+k] != NULL) {
//		renderTexture(num_tex[j*npressure+k], renderer, x + j*d, y - k*d - 13);
//	}

	SDL_RenderPresent(renderer);
	SDL_Delay(16);

	return 0;
}

//-------------------------------------------------------------------
//                      Click handling subroutine
//-------------------------------------------------------------------

//int handleClick(SDL_Event e, int ntemp,
//		SDL_Surface **pies, int *xstart, int *xend, int *ystart, int *yend,
//		SDL_Texture **background_tex, double *pie_radius, char *path, SDL_Texture ***num_tex) {
//
//	int i = 0;
//	int x = 0; int y = 0;
//	int xWR = 0; int yWR = 0;
//	int xtopic = 0; int ytopic = 0;
//	Uint32 *pixmem32;
//
//	// Reset screen
//	for (x=0;x<(*pies)->w;x++) {
//		for (y=0;y<=(*pies)->h;y++) {
//			pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
//			*pixmem32 = SDL_MapRGBA((*pies)->format, 0, 0, 0, 0);
//		}
//	}
//	for (i=0;i<nnum;i++) (*num_tex)[i] = NULL;
//
//	// Dummy response to click
//	if (e.button.x >= 4 && e.button.x <= 64 && e.button.y >= 444 && e.button.y <= 491) i = 0;          // Potassium
//
//	return 0;
//}

//-------------------------------------------------------------------
//                    Number plotting subroutine
//-------------------------------------------------------------------

//int PlotNumChem(int PT, int ntemp, double Tmin, double Tstep, int npressure, double Pmin, double Pstep, int npH,
//		double pHmin, double pHstep, int npe, double pemin, double pestep, int nWR, double WRmin, double WRstep,
//		SDL_Texture ***Numbers, char* FontFile) {
//
//	int i = 0;
//	char nb[20];
//	SDL_Color black;
//	black.r = 0; black.g = 0; black.b = 0; black.a = 0;
//	SDL_Color white;
//	white.r = 255; white.g = 255; white.b = 255; white.a = 0;
//
//	for (i=0;i<ntemp;i++) {
//		if (i == 0 && Tmin == 0) scanNumber(&nb, 0.01);    // Right-justified
//		else scanNumber(&nb, Tmin + (double) i*Tstep);     // Right-justified
//		(*Numbers)[i] = renderText(nb, FontFile, black, 14, renderer);
//	}
//
//	return 0;
//}

//-------------------------------------------------------------------
//                     Chamber plot subroutine
//-------------------------------------------------------------------

int ChamberPlot(SDL_Surface **chamber, char *FontFile, int ntemp, double radius, double **R2Hcompo, double Rchamber, double R1, int nsalts, int contour, int dbase_type) {

	SDL_Color white;
	SDL_Color aqua;
	SDL_Color cyan;
	white.r = 250; white.g = 250; white.b = 250; white.a = 255;
	aqua.r = 0; aqua.g = 128; aqua.b = 255; aqua.a = 255;
	cyan.r = 138; cyan.g = 240; cyan.b = 255; cyan.a = 255;

	SDL_Color color;

	int x = 0; int y = 0; // Chamber center coordinates
	int xvar = 0; int yvar = 0;
	int i = 0; int j = 0; int k = 0; int l = 0;
	Uint32 *pixmem32;

	int r = 0; int g = 0; int b = 0; int a = 0;

	x = 400; y = 325;

	// Outer chamber
	for (i=0;i<2*(int)radius;i++) {
		for (j=0;j<2*(int)radius;j++) {

			xvar = x - (int)radius + i; yvar = y - (int)radius + j;

			color = cyan;
			r = color.r; g = color.g; b = color.b; a = ((int)R2Hcompo[1][0]-250)/23.15*color.a;
			if (sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) < radius) {
				pixmem32 = (Uint32*) (*chamber)->pixels + yvar*(*chamber)->w + xvar;
				*pixmem32 = SDL_MapRGBA((*chamber)->format, (r*(1-abs(y-yvar)/radius) + 2*r)/3,
														    (g*(1-abs(y-yvar)/radius) + 2*g)/3,
														    (b*(1-abs(y-yvar)/radius) + 2*b)/3, a);
			}
			color = aqua;
			r = color.r; g = color.g; b = color.b; a = color.a;
			if (sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) < radius*R1/Rchamber) {
				pixmem32 = (Uint32*) (*chamber)->pixels + yvar*(*chamber)->w + xvar;
				*pixmem32 = SDL_MapRGBA((*chamber)->format, (r*(1-abs(y-yvar)/radius) + 2*r)/3,
														    (g*(1-abs(y-yvar)/radius) + 2*g)/3,
														    (b*(1-abs(y-yvar)/radius) + 2*b)/3, a);
			}
		}
	}
	// Inner salt layers
	for (k=1;k<ntemp;k++) {
		if (R2Hcompo[k][0] <= 0.0) break;
		for (i=0;i<2*(int)radius;i++) {
			for (j=0;j<2*(int)radius;j++) {

				xvar = x - (int)radius + i; yvar = y - (int)radius + j;

				// Ice rim
				color = cyan;
				r = color.r; g = color.g; b = color.b; a = ((int)R2Hcompo[k][0]-250)/23.15*color.a;
				if (sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) < radius*R2Hcompo[k-1][1]*R1/Rchamber
						&& sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) >= radius*R2Hcompo[k][1]*R1/Rchamber
						&& yvar-y <= (int)(radius*R1/Rchamber*(1.0-R2Hcompo[k][2]))) {
					pixmem32 = (Uint32*) (*chamber)->pixels + yvar*(*chamber)->w + xvar;
					*pixmem32 = SDL_MapRGBA((*chamber)->format, (r*(1-abs(y-yvar)/radius) + 2*r)/3,
															    (g*(1-abs(y-yvar)/radius) + 2*g)/3,
															    (b*(1-abs(y-yvar)/radius) + 2*b)/3, a);
				}
				// Salt deposit
				double gradient = (double)(xvar-(x-ceil(radius*R2Hcompo[k][1]*R1/Rchamber)))
						                         / (2.0*radius*R2Hcompo[k][1]*R1/Rchamber);
				double pos = 0.0;
				double pos_old = 0.0;
				for (l=0;l<nsalts;l++) {
					pos_old = pos;
					pos += R2Hcompo[k][l+3];
					if (gradient <= pos/100.0 && gradient > pos_old/100.0) color = minColor(R2Hcompo[0][l+3], dbase_type);
				}

				r = color.r; g = color.g; b = color.b; a = color.a;
				if (sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) < radius*R2Hcompo[k-1][1]*R1/Rchamber
						&& yvar-y > (int)(radius*R1/Rchamber*(1.0-R2Hcompo[k][2]))
						&& yvar-y <= (int)(radius*R1/Rchamber*(1.0-R2Hcompo[k-1][2]))) {
					pixmem32 = (Uint32*) (*chamber)->pixels + yvar*(*chamber)->w + xvar;
					*pixmem32 = SDL_MapRGBA((*chamber)->format, (r*(1-abs(y-yvar)/radius) + 2*r)/3,
															    (g*(1-abs(y-yvar)/radius) + 2*g)/3,
															    (b*(1-abs(y-yvar)/radius) + 2*b)/3, a);
				}
				// Salt deposit contour
				if (contour) {
					color = white;
					r = color.r; g = color.g; b = color.b; a = color.a;
					if (((int)(sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y))) == (int)(radius*R2Hcompo[k-1][1]*R1/Rchamber)
							&& yvar-y >= (int)(radius*R1/Rchamber*(1.0-R2Hcompo[k][2]))
							&& yvar-y <= (int)(radius*R1/Rchamber*(1.0-R2Hcompo[k-1][2])))
					 || (sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) <= radius*R2Hcompo[k-1][1]*R1/Rchamber && yvar-y == (int)(radius*R1/Rchamber*(1.0-R2Hcompo[k-1][2])))) {
						pixmem32 = (Uint32*) (*chamber)->pixels + yvar*(*chamber)->w + xvar;
						*pixmem32 = SDL_MapRGBA((*chamber)->format, (r*(1-abs(y-yvar)/radius) + 2*r)/3,
																	(g*(1-abs(y-yvar)/radius) + 2*g)/3,
																	(b*(1-abs(y-yvar)/radius) + 2*b)/3, a);
					}
				}
			}
		}
	    // If ice covers salt, reset salt deposit height
	    if (R2Hcompo[k][1] < 1.0-R2Hcompo[k][2]) {
	    	printf("R2/R1=%g < 1-H/R1=1-%g, resetting R1 and H at k=%d, T=%g\n", R2Hcompo[k][1], R2Hcompo[k][2], k, R2Hcompo[k][0]);
	        R1 *= R2Hcompo[k][1];
	        double R2mem = R2Hcompo[k][1];
	        double Hmem = R2Hcompo[k][2];
	        for (i=0;i<ntemp;i++) {
	        	R2Hcompo[i][1] /= R2mem;
	        	R2Hcompo[i][2] -= Hmem;
	        	R2Hcompo[i][2] /= R2mem;
	        }
//	        printf("T(K) \tR2/R1 \tH/R1\n");
//	        for (i=k;i<ntemp;i++) printf("%g \t %g \t %g\n", R2Hcompo[i][0], R2Hcompo[i][1], R2Hcompo[i][2]);
	    }
	}

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine ExtractWrite
 *
 * Write selected output from PHREEQC
 *
 *--------------------------------------------------------------------*/

int ExtractWrite(int instance, double*** data, int line, int nvar) {
	VAR v;
	int i = 0;
	VarInit(&v);

	for (i=0;i<nvar;i++) {
		GetSelectedOutputValue(instance,1,i,&v);
		if (fabs(v.dVal) < 1e-50) (*data)[line][i] = 0.0;
		else (*data)[line][i] = v.dVal;
	}
	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine MineralName
 *
 * Name mineral based on index
 *
 *--------------------------------------------------------------------*/

int MineralName(int i, char **name, int dbase_type) {

	strcpy(*name, "Name not found, check indices in MineralName()");
	int j = 0;
	if (!dbase_type) j = 15; // TC17 database
	else j = 18;             // Carbonate-Si-CH4 database
	if (i == j) strcpy(*name, "Ice(s)"); j = j+2;
	if (i == j) strcpy(*name, "Halite"); j = j+2;
	if (i == j) strcpy(*name, "Hydrohalite"); j = j+2;
	if (i == j) strcpy(*name, "Sylvite"); j = j+2;
	if (i == j) strcpy(*name, "Antarcticite"); j = j+2;
	if (i == j) strcpy(*name, "Bischofite"); j = j+2;
	if (i == j) strcpy(*name, "MgCl2:8H2O"); j = j+2;
	if (i == j) strcpy(*name, "MgCl2:12H2O"); j = j+2;
	if (i == j) strcpy(*name, "Carnallite"); j = j+2;
	if (i == j) strcpy(*name, "Tachyhydrite"); j = j+2;
	if (i == j) strcpy(*name, "Anhydrite"); j = j+2;
	if (i == j) strcpy(*name, "Arcanite"); j = j+2;
	if (i == j) strcpy(*name, "Epsomite"); j = j+2;
	if (i == j) strcpy(*name, "Gypsum"); j = j+2;
	if (i == j) strcpy(*name, "Meridianiite"); j = j+2;
	if (i == j) strcpy(*name, "Mirabilite"); j = j+2;
	if (i == j) strcpy(*name, "Picromerite"); j = j+2;
	if (i == j) strcpy(*name, "Glaserite"); j = j+2;
	if (i == j) strcpy(*name, "Bloedite"); j = j+2;
	if (i == j) strcpy(*name, "Kieserite"); j = j+2;
	if (i == j) strcpy(*name, "Thenardite"); j = j+2;
	if (i == j) strcpy(*name, "Hexahydrite"); j = j+2;

	if (!dbase_type) { // TC17 database
		if (i == j) strcpy(*name, "Na2SO4:7H2O"); j = j+2;
		if (i == j) strcpy(*name, "Starkeyite"); j = j+2;
		if (i == j) strcpy(*name, "Pentahydrite"); j = j+2;
		if (i == j) strcpy(*name, "Bassanite"); j = j+2;
		if (i == j) strcpy(*name, "Glauberite"); j = j+2;
		if (i == j) strcpy(*name, "Labile_Salt"); j = j+2;
		if (i == j) strcpy(*name, "Leonite"); j = j+2;
		if (i == j) strcpy(*name, "Langbeinite"); j = j+2;
		if (i == j) strcpy(*name, "Syngenite"); j = j+2;
		if (i == j) strcpy(*name, "Polyhalite"); j = j+2;
		if (i == j) strcpy(*name, "Kainite");
	}
	else { // carbonate-Si-CH4 database
		if (i == j) strcpy(*name, "Aragonite"); j = j+2;
		if (i == j) strcpy(*name, "Artinite"); j = j+2;
		if (i == j) strcpy(*name, "Calcite"); j = j+2;
		if (i == j) strcpy(*name, "Dolomite"); j = j+2;
		if (i == j) strcpy(*name, "Hydromagnesite"); j = j+2;
		if (i == j) strcpy(*name, "Ikaite"); j = j+2;
		if (i == j) strcpy(*name, "Kalicinite"); j = j+2;
		if (i == j) strcpy(*name, "Landsfordite"); j = j+2;
		if (i == j) strcpy(*name, "Magnesite"); j = j+2;
		if (i == j) strcpy(*name, "Na2CO3:7H2O"); j = j+2;
		if (i == j) strcpy(*name, "Nahcolite"); j = j+2;
		if (i == j) strcpy(*name, "Natron"); j = j+2;
		if (i == j) strcpy(*name, "Nesquehonite"); j = j+2;
		if (i == j) strcpy(*name, "Trona"); j = j+2;
		if (i == j) strcpy(*name, "Vaterite"); j = j+2;
		if (i == j) strcpy(*name, "Akermanite"); j = j+2;
		if (i == j) strcpy(*name, "Anthophyllite"); j = j+2;
		if (i == j) strcpy(*name, "Antigorite"); j = j+2;
		if (i == j) strcpy(*name, "Chalcedony"); j = j+2;
		if (i == j) strcpy(*name, "Chrysotile"); j = j+2;
		if (i == j) strcpy(*name, "Diopside"); j = j+2;
		if (i == j) strcpy(*name, "Enstatite"); j = j+2;
		if (i == j) strcpy(*name, "Forsterite"); j = j+2;
		if (i == j) strcpy(*name, "Quartz"); j = j+2;
		if (i == j) strcpy(*name, "Sepiolite"); j = j+2;
		if (i == j) strcpy(*name, "SiO2(am)"); j = j+2;
		if (i == j) strcpy(*name, "Talc"); j = j+2;
		if (i == j) strcpy(*name, "Methane_hydrate");
	}
	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine Vm
 *
 * Return mineral volume based on name by reading PHREEQC database
 *
 *--------------------------------------------------------------------*/

double Vm(char *name, char *dbase) {

	double molarvol = 0.0;

	int i = 0;
	int same = 0;
	FILE *f;
	int line_length = 1024;
	char line[line_length]; // Individual line
	line[0] = '\0';

	f = fopen (dbase,"r"); 	// Open PHREEQC database file
	if (f == NULL) printf("Vm: Missing database file. Path: %s\n", dbase);

	while (fgets(line, line_length, f)) {
		for (i=0;i<20;i++) {
			if (i == 4) {
				if (line[i] != '2' && line[i] != 'o') {
					same = 1;
					break;
				}
			}
			if (i == 6) {
				same = 1;
				break; // To distinguish hydrohalite and hydromagnesite, of MgCl2:8H2O and :12H2O
			}
			if (line[i] != name[i]) break;
		}
		if (same && line[1] == '-' && line[2] == 'a' && line[3] == 'n' && line[4] == 'a' ) { // -analytic line has been read, -Vm is upcoming line
			fseek(f, 5, SEEK_CUR);
			fscanf(f, "%lg", &molarvol);
			break;
		}
	}
	return molarvol;
}

/*--------------------------------------------------------------------
 *
 * Subroutine VsphSegm
 *
 * Return the volume of a spherical segment of radius R2, height x,
 * and whose far side is at distance b from center of sphere
 *
 *--------------------------------------------------------------------*/

double VsphSegm(double x, double R2, double b) {

	double segm = 0.0;

	segm = (R2*R2 - b*b)*x + b*x*x - x*x*x/3.0;
    segm *= M_PI;

	return segm;
}

/*--------------------------------------------------------------------
 *
 * Subroutine R2iceCap
 *
 * Return the volume of a hollow spherical cap encasing volume
 * Vsol, of inner radius x and outer radius R1, and truncated by height H
 *
 *--------------------------------------------------------------------*/

double R2iceCap (double x, double Vsol, double R1, double H) {

	double FindR2bigH = 0.0;
    FindR2bigH = Vsol/M_PI - 2.0/3.0*x*x*x + x*x*(H-R1) + R1*R1*R1/3.0 - H*H*H/3.0 - R1*R1*H + H*H*R1;

	return FindR2bigH;
}

/*--------------------------------------------------------------------
 *
 * Subroutine minColor
 *
 * Color-codes minerals on chamber plot
 *
 *--------------------------------------------------------------------*/

SDL_Color minColor (int i, int dbase_type) {

	SDL_Color color;

	SDL_Color cyan;			cyan.r = 138; cyan.g = 240; cyan.b = 255; cyan.a = 255;

	SDL_Color orchid;		orchid.r = 122; orchid.g = 129; orchid.b = 255; orchid.a = 255;
	SDL_Color lavender;		lavender.r = 215; lavender.g = 131; lavender.b = 255; lavender.a = 255;
	SDL_Color grape;		grape.r = 148; grape.g = 55; grape.b = 255; grape.a = 255;
	SDL_Color purple;		purple.r = 168; purple.g = 50; purple.b = 208; purple.a = 255;
	SDL_Color eggplant;		eggplant.r = 83; eggplant.g = 27; eggplant.b = 147; eggplant.a = 255;
	SDL_Color plum;			plum.r = 148; plum.g = 33; plum.b = 147; plum.a = 255;

	SDL_Color bubblegum;	bubblegum.r = 255; bubblegum.g = 133; bubblegum.b = 255; bubblegum.a = 255;
	SDL_Color carnation;	carnation.r = 255; carnation.g = 138; carnation.b = 216; carnation.a = 255;
	SDL_Color magenta;		magenta.r = 255; magenta.g = 64; magenta.b = 255; magenta.a = 255;
	SDL_Color pink;			pink.r = 255; pink.g = 47; pink.b = 146; pink.a = 255;

	SDL_Color salmon;		salmon.r = 255; salmon.g = 126; salmon.b = 121; salmon.a = 255;
	SDL_Color maraschino;	maraschino.r = 255; maraschino.g = 38; maraschino.b = 0; maraschino.a = 255;
	SDL_Color red;			red.r = 250; red.g = 20; red.b = 20; red.a = 255;
	SDL_Color cayenne;		cayenne.r = 148; cayenne.g = 17; cayenne.b = 0; cayenne.a = 255;
	SDL_Color maroon;		maroon.r = 128; maroon.g = 0; maroon.b = 64; maroon.a = 255;
	SDL_Color maroon2;		maroon2.r = 148; maroon2.g = 23; maroon2.b = 81; maroon2.a = 255;

	SDL_Color cantaloupe;	cantaloupe.r = 255; cantaloupe.g = 212; cantaloupe.b = 121; cantaloupe.a = 255;
	SDL_Color tangerine;	tangerine.r = 255; tangerine.g = 147; tangerine.b = 0; tangerine.a = 255;
	SDL_Color orange;		orange.r = 238; orange.g = 124; orange.b = 22; orange.a = 255;
	SDL_Color mocha;		mocha.r = 148; mocha.g = 82; mocha.b = 0; mocha.a = 255;

	SDL_Color banana;		banana.r = 255; banana.g = 252; banana.b = 121; banana.a = 255;
	SDL_Color lemon;		lemon.r = 255; lemon.g = 251; lemon.b = 0; lemon.a = 255;
	SDL_Color yellow;		yellow.r = 245; yellow.g = 217; yellow.b = 33; yellow.a = 255;
	SDL_Color gold;			gold.r = 255; gold.g = 255; gold.b = 158; gold.a = 255;
	SDL_Color asparagus;	asparagus.r = 146; asparagus.g = 144; asparagus.b = 0; asparagus.a = 255;

	SDL_Color ice;			ice.r = 115; ice.g = 253; ice.b = 255; ice.a = 255;
	SDL_Color sky;			sky.r = 118; sky.g = 214; sky.b = 255; sky.a = 255;
	SDL_Color sea_foam;		sea_foam.r = 0; sea_foam.g = 250; sea_foam.b = 146; sea_foam.a = 255;
	SDL_Color turquoise;	turquoise.r = 0; turquoise.g = 253; turquoise.b = 255; turquoise.a = 255;

	SDL_Color aqua;			aqua.r = 0; aqua.g = 128; aqua.b = 255; aqua.a = 255;
	SDL_Color blueberry;	blueberry.r = 4; blueberry.g = 51; blueberry.b = 255; blueberry.a = 255;
	SDL_Color blue;			blue.r = 0; blue.g = 0; blue.b = 255; blue.a = 255;
	SDL_Color ocean;		ocean.r = 0; ocean.g = 84; ocean.b = 147; ocean.a = 255;
	SDL_Color midnight;		midnight.r = 1; midnight.g = 25; midnight.b = 147; midnight.a = 255;
	SDL_Color teal;			teal.r = 0; teal.g = 145; teal.b = 147; teal.a = 255;

	SDL_Color honeydew;		honeydew.r = 212; honeydew.g = 251; honeydew.b = 121; honeydew.a = 255;
	SDL_Color lime;			lime.r = 142; lime.g = 250; lime.b = 0; lime.a = 255;
	SDL_Color fern;			fern.r = 79; fern.g = 143; fern.b = 0; fern.a = 255;
	SDL_Color flora;		flora.r = 115; flora.g = 250; flora.b = 121; flora.a = 255;
	SDL_Color spring;		spring.r = 0; spring.g = 249; spring.b = 0; spring.a = 255;
	SDL_Color clover;		clover.r = 0; clover.g = 143; clover.b = 0; clover.a = 255;
	SDL_Color green;		green.r = 39; green.g = 145; green.b = 39; green.a = 255;
	SDL_Color light_green;	light_green.r = 204; light_green.g = 255; light_green.b = 102; light_green.a = 255;
	SDL_Color spindrift;	spindrift.r = 102; spindrift.g = 255; spindrift.b = 204; spindrift.a = 255;
	SDL_Color moss;			moss.r = 0; moss.g = 143; moss.b = 0; moss.a = 255;

	SDL_Color white; 		white.r = 250; white.g = 250; white.b = 250; white.a = 255;
	SDL_Color gray;			gray.r = 174; gray.g = 174; gray.b = 174; gray.a = 255;
	SDL_Color iron;			iron.r = 94; iron.g = 94; iron.b = 94; iron.a = 255;
	SDL_Color black; 		black.r = 30; black.g = 30; black.b = 30; black.a = 255;

	int j = 0;
	if (!dbase_type) j = 15; // TC17 database
	else j = 18; // Carbonate-Si-CH4 database

	if (i == j) color = cyan;
	j = j+2;
	if (i == j) color = orchid;
	j = j+2;
    if (i == j) color = lavender;
    j = j+2;
    if (i == j) color = grape;
    j = j+2;
    if (i == j) color = purple;
    j = j+2;
    if (i == j) color = eggplant;
    j = j+2;
    if (i == j) color = plum;
    j = j+2;
    if (i == j) color = bubblegum;
    j = j+2;
    if (i == j) color = carnation;
    j = j+2;
    if (i == j) color = magenta;
    j = j+2;
    if (i == j) color = pink;
    j = j+2;
    if (i == j) color = banana;
    j = j+2;
    if (i == j) color = lemon;
    j = j+2;
    if (i == j) color = yellow;
    j = j+2;
    if (i == j) color = gold;
    j = j+2;
	if (i == j) color = cantaloupe;
	j = j+2;
	if (i == j) color = tangerine;
	j = j+2;
    if (i == j) color = orange;
    j = j+2;
    if (i == j) color = mocha;
    j = j+2;
    if (i == j) color = salmon;
    j = j+2;
    if (i == j) color = maraschino;
    j = j+2;
    if (i == j) color = red;
    j = j+2;
	if (i == j) color = blue;
	j = j+2;
	if (i == j) color = aqua;
	j = j+2;
	if (i == j) color = turquoise;
	j = j+2;
	if (i == j) color = ocean;
	j = j+2;
    if (i == j) color = cayenne;
    j = j+2;
    if (i == j) color = teal;
    j = j+2;
	if (i == j) color = blueberry;
	j = j+2;
    if (i == j) color = maroon;
    j = j+2;
    if (i == j) color = maroon2;
    j = j+2;
    if (i == j) color = ice;
    j = j+2;
    if (i == j) color = sea_foam;
    j = j+2;
    if (i == j) color = sky;
    j = j+2;
	if (i == j) color = midnight;
	j = j+2;
    if (i == j) color = honeydew;
    j = j+2;
    if (i == j) color = lime;
    j = j+2;
    if (i == j) color = fern;
    j = j+2;
    if (i == j) color = flora;
    j = j+2;
    if (i == j) color = spring;
    j = j+2;
    if (i == j) color = clover;
    j = j+2;
    if (i == j) color = green;
    j = j+2;
    if (i == j) color = light_green;
    j = j+2;
    if (i == j) color = spindrift;
    j = j+2;
    if (i == j) color = moss;
    j = j+2;
    if (i == j) color = asparagus;
    j = j+2;
    if (i == j) color = white;
    j = j+2;
    if (i == j) color = gray;
    j = j+2;
    if (i == j) color = iron;
    j = j+2;
    if (i == j) color = black;
    j = j+2;

	return color;
}
