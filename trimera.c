
//TRIMER HELIJA Josip
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c"

#define Ns 500     // broj koraka
#define Nw 500     // broj setaca
#define Nb 10    // broj blokova
#define NbSkip 1 // broj prvih blokova koje preskacemo




#define R0 7.537   // Angstorm, širina jame
#define V0 565.44 // mK, dubina jame
#define D 6059.64 //D je konstanta h^2/(2*m*amu*k_b)
//https://www.wolframalpha.com/input?i=+%286.582119569e-16%29%5E2%2F%282*4.002602*%28931.49410242e6%2F%28299792458*1e10%29%5E2%29*8.617333262e-8%29
// kad se uvrste sve mjerne jedinice i pokrate, energiju dobijemo u mK.


// Varijacijski parametri RoVo jamu:
double kij, Rij;
void fprovjera(FILE *dat)
{
    if (dat == NULL)
    {
        printf("Error while opening the file.\n");
        exit(EXIT_FAILURE);
    }
}

double Psi(double r)//literatura sa one drive (4.31).
{
    // RoVo
     if(r <= Rij){
        return sin(kij*r)/r;
     }
     else{
        return (sin(kij*Rij)/r)*exp((kij*(r-Rij))/tan(kij*Rij));
     }
}

double fdr(double r)//literatura 4.97a)
{
    if(r <= Rij){
        return (kij*r*pow(tan(kij*r),-1)-1.0)/pow(r,2);
    }
    else{
        return (kij*r*pow(tan(kij*Rij),-1)-1.0)/pow(r,2);
    }
}

double fddr(double r)// 4.97b)
{
    if(r <= Rij ){
        return (kij*r*(2.0*pow(tan(kij*r),-1)-kij*r*pow(sin(kij*r),-2))- 1.0)/pow(r,2);
    }
    else{
        return (2.0*kij*r*pow(tan(kij*Rij),-1)-1.0)/pow(r,2);
    }
}


double energija(double r[4][4], double x[4][4][4])//local energy (4.96)
{
    // RoVo potencijal za trimer
     double V12=0,V13=0,V23=0,Vuk=0;

     if(r[1][2] <= R0){
        V12 = -V0;
        }
     if(r[1][3] <= R0){
        V13 = -V0;
        }
     if(r[2][3] <= R0){
        V23=-V0;
        }
Vuk=V12+V13+V23;
     return  Vuk-D*(fddr(r[1][2])+fddr(r[1][3])+(pow((fdr(r[1][2])*x[1][2][1]+fdr(r[1][3])*x[1][3][1]),2)+pow((fdr(r[1][2])*x[1][2][2]+fdr(r[1][3])*x[1][3][2]),2)+pow((fdr(r[1][2])*x[1][2][3]+fdr(r[1][3])*x[1][3][3]),2))+fddr(r[2][1])+fddr(r[2][3])+(pow((fdr(r[2][1])*x[2][1][1]+fdr(r[2][3])*x[2][3][1]),2)+pow((fdr(r[2][1])*x[2][1][2]+fdr(r[2][3])*x[2][3][2]),2)+pow((fdr(r[2][1])*x[2][1][3]+fdr(r[2][3])*x[2][3][3]),2))+fddr(r[3][1])+fddr(r[3][2])+(pow((fdr(r[3][1])*x[3][1][1]+fdr(r[3][2])*x[3][2][1]),2)+pow((fdr(r[3][1])*x[3][1][2]+fdr(r[3][2])*x[3][2][2]),2)+pow((fdr(r[3][1])*x[3][1][3]+fdr(r[3][2])*x[3][2][3]),2)));
   //return V12+V13+V23;
}

#include "VMC.c"

int main()
{
     //double Rrange[] = {9.1,10.5}, dR;             // Rij (6.95) i (0.23,0.28).minimum
     Rij=9.5;
     kij=0.185;                      // 14.6;
     //double k = sqrt((V0*(1./6059.64))); // valni vektor iz potencijala V0 za pravokutnu jamu k=0.305
     //double krange[] = {0.165, 0.23}, dk;            // kij ,0.229,0.259  0.16, 0.20} 10.1;


    double E[2];      // Energija[1]  i greška[2]
    double min_par[2][4]; // 1. indeks prva i druga najmanja E. 2.indeks je: Rij-[0], kij-[1], E-[2], sigmaE-[3].
    //int N =9;        // resetka parametara NxN



    FILE *param;
    param = fopen("svi_parametri9.5.txt", "w");
    fprovjera(param);
/*
    // koraci
    dR = (Rrange[1] - Rrange[0]) / N;
    dk = (krange[1] - krange[0]) / N;

    int ir, ik, i, j;

    for (i = 0; i < 2; i++){
        for (j = 0; j < 4; j++){
            min_par[i][j] = 0.0;
            }
    }

    for (ir = 0; ir <= N; ir++)
    {
        Rij = Rrange[0] + ir * dR;

    for (ik = 0; ik <= N; ik++)
        {
            kij = krange[0] + ik * dk;*/

            VMC(E);

            fprintf(param, "%e\t%e\t%e\t%e\n", Rij, kij, E[0], E[1]);


        /*if (E[0] < min_par[0][2]) // manji od najmanje energije
            {

                for (i = 0; i < 4; i++) //  pomak 1. na 2. mjesto
                {
                    min_par[1][i] = min_par[0][i];
                }*/
                min_par[0][0] = Rij;
                min_par[0][1] = kij;
                min_par[0][2] = E[0];//energija
                min_par[0][3] = E[1];//greska
           // }


        //}
        fprintf(param, "\n\n");

    //}

fclose(param);

    param = fopen("minimalni_paramteri9.5.txt", "w");
    fprovjera(param);
    fprintf(param, "Najmanje energije sa parametrima:\n");

        fprintf(param, " Rij = %lf, kij = %lf\t E = %8.5e +- %6.2e\n",  min_par[0][0], min_par[0][1], min_par[0][2], min_par[0][3]);

    fclose(param);

    return 0;
}
