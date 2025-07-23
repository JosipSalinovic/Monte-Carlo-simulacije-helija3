


//TRIMER HELIJA Josip
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c"

const long int  Ns=1000;     // broj koraka
const long int  Nw= 500;     // broj setaca
const long int  Nb= 10;   // broj blokova
const long int  NbSkip= 1; // broj prvih blokova koje preskacemo





const double  Dk= 6059.6499694573; //D je konstanta h^2/(2*m*amu*k_b)
//https://www.wolframalpha.com/input?i=+%286.582119569e-16%29%5E2%2F%282*4.002602*%28931.49410242e6%2F%28299792458*1e10%29%5E2%29*8.617333262e-8%29
// kad se uvrste sve mjerne jedinice i pokrate, energiju dobijemo u mK.
// HFDB konstante
const double  A= 184431.01;
const double  epsilon= 10948; // mK
const double  rm= 2.963;      // Å
const double  alfa= 10.43329537;
const double  beta= -2.27965105;
const double  C6= 1.36745214;
const double  C8= 0.42123807;
const double  C10= 0.17473318;
const double  D= 1.4826;


double g;
double s;
double a;
void fprovjera(FILE *dat)
{
    if (dat == NULL)
    {
        printf("Error while opening the file.\n");
        exit(EXIT_FAILURE);
    }
}
double F(double x)
{
    if(x<=D){
        return exp(-pow(((D/x)-1.),2));}
    else{
        return 1.0;}
}
double V(double r){

    return epsilon*(A*exp(-alfa*(r/rm)+beta*pow(r/rm,2))-F(r/rm)*(C6/pow(r/rm,6) + C8/pow(r/rm,8)+C10/pow(r/rm,10)));



}
double Psi(double r)//f2
{
    //return exp(-g*exp(-a*r)-s*r);

     return (exp(-pow((a/r),g)-s*r))/r;
}

double fdr(double r)//literatura 4.91
{
        //return (a*g*exp(-a*r)-s)/r;
        return (-s*r+g*pow((a/r),g)-1.0)/pow(r,2);


}

double fddr(double r)// 4.94
{
     //return (-2*s+(2*a*g-a*a*g*r)*exp(-a*r))/r;
     return (-2.0*s*r-1.000+(-pow(g,2)+g)*pow((a/r),g))/pow(r,2);

}


double energija(double r[4][4], double x[4][4][4])//local energy (4.96)
{
    double Vuk=(V(r[1][2])+V(r[1][3])+V(r[2][3]));

     return  Vuk-Dk*(fddr(r[1][2])+fddr(r[1][3])+(pow((fdr(r[1][2])*x[1][2][1]+fdr(r[1][3])*x[1][3][1]),2)+pow((fdr(r[1][2])*x[1][2][2]+fdr(r[1][3])*x[1][3][2]),2)+pow((fdr(r[1][2])*x[1][2][3]+fdr(r[1][3])*x[1][3][3]),2))+fddr(r[2][1])+fddr(r[2][3])+(pow((fdr(r[2][1])*x[2][1][1]+fdr(r[2][3])*x[2][3][1]),2)+pow((fdr(r[2][1])*x[2][1][2]+fdr(r[2][3])*x[2][3][2]),2)+pow((fdr(r[2][1])*x[2][1][3]+fdr(r[2][3])*x[2][3][3]),2))+fddr(r[3][1])+fddr(r[3][2])+(pow((fdr(r[3][1])*x[3][1][1]+fdr(r[3][2])*x[3][2][1]),2)+pow((fdr(r[3][1])*x[3][1][2]+fdr(r[3][2])*x[3][2][2]),2)+pow((fdr(r[3][1])*x[3][1][3]+fdr(r[3][2])*x[3][2][3]),2)));
}


#include "VMChfdb2.c"
int main(){


 double grange[] = {4.165, 4.17}, dg; // g
    double arange[] = {2.815,2.81}, da; // a
     double srange[] = {0.0245,0.0250}, ds; // s

    double E[2];      // E i sigma E

    int N = 1;        // N^2 točaka u parametarskom prostoru se provjerava



    FILE *param;
    param = fopen("parametri2.txt", "w");
    fprovjera(param);

    // koraci u parametarskom prostoru
    dg = (grange[1] - grange[0]) / N;
    da = (arange[1] - arange[0]) / N;
    ds = (srange[1] - srange[0]) / N;

    int ig, ia, i, j,is;

    for (is = 0; is < N; is++){
            s=srange[0]+is*ds;

         for (ig = 0; ig < N; ig++){

             g = grange[0] + ig * dg;


              for (ia = 0; ia < N; ia++){
                    a = arange[0] + ia * da;

                    VMC(E);
                   fprintf(param, "%e\t%e\t%e\t%e\t%e\n", g, a, s,E[0],E[1]);
        }}}






fclose(param);
}


