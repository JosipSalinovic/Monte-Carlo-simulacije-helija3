/* ====================================================================
   _______________________ ______  __   ______________________  ___
   __  ___/___  __/__  __ \___  / / /   __  ___/____  _/___   |/  /
   _____ \ __  /   _  / / /__  /_/ /    _____ \  __  /  __  /|_/ /
   ____/ / _  /    / /_/ / _  __  /     ____/ / __/ /   _  /  / /
   /____/  /_/     \____/  /_/ /_/      /____/  /___/   /_/  /_/
	  Prirodoslovno-matematicki fakultet u Splitu
	  Stohasticke simulacije u klasicnoj i kvantnoj fizici
	  Varijacijski Monte Carlo :: Osnovno stanje vodikova atoma
	  2018/2019, 2019/2020
Koristene oznake:
	  Nw    = broj setaca (number of walkers)
	  iw    = indeks setaca (index of walker)
	  Nb    = broj blokova (number of blocks)
	  ib    = indeks bloka (index of block)
	  Nt    = broj koraka (number of time-steps)
	  it    = indeks vremenskog koraka (time-step)
	  k     = koordinata (coordinate)
	  NbSkip= broj preskocenih blokova
	  E     = lokalna energija (local energy)
	  SwE   = suma (srednjih E) po setacima
	  StE   = suma (srednjih E) po koracima
	  SbE   = suma (srednjih E) po blokovima
	  accept= brojac prihvacenih koraka
	  acc_ib= udio prihvacenih koraka
	  x     = koordinate posa setaca
	  dx    = promjena koordinate nasumicno od -dk do dk
	  dk    = maksimalna duljina koraka
	  xp    = koordinate probnog posa setaca
	  r1    = modul radijvektora probnog posa
	  r2    = r1*r1
	  P     = vjerojatnost nalazenja na posu x
	  Pp    = vjerojatnost nalazenja u probnom posu
	  T     = vjerojatnost prijelaza x -> xp
	  Psi   = valna funkcija
======================================================================= */

void VMC(double *E_ret){
	long idum = (-58755);
	int is, ib, iw, k, Nbskip=1, itmp,a,d;
	double accept, acc_is, AE, sigmaE;
	double x[4][4][4], xp[4][4][4], dk[4], E[Nw + 1], P[Nw + 1], pos[4][4][Nw],newpos[4][4];
	double dx, r[4][4],rnew[4][4], SwE, SsE, SbE, SbE2, Pp, T;

	accept = 0.;				 // prihvacanje
	dk[1] = dk[2] = dk[3] = 3.2; // maksimalne promjene koordinata
	FILE *fout;


		char ime[100];
		 sprintf(ime, "E_RoVo_R_%.3lf_k_%.3lf.dat", Rij, kij);
		//sprintf(ime, "E_HFDB_s_%lf_b_%lf.dat", s, b);
		fout = fopen(ime, "w"); // datoteka za pohranu srednjih vrijednosti
		fprovjera(fout);


	// inicijalizacija posa gdje je gustoca Psi*Psi znacajna
	for (iw = 1; iw <= Nw; iw++)
	{
		r[1][2] =0.;
		r[1][3] =0.;
		r[2][3] =0.;
		r[2][1] =0.;
		r[3][1] =0.;
		r[3][2] =0.;


	 for(a=1;a<=3;a++){// petlje po atomima

			 for(k=1;k<=3;k++){

                  pos[a][k][iw]= 60.0 * (ran1(&idum) - 0.5);// koordinate za tri atoma
			 }}

			 for(k=1;k<=3;k++){

                //dvočestične komponenete. x[a][d][k] gdje je [a] = prvi atom [d]= drugi atom, [k]-komponenta x,y,z.

                x[1][2][k]=pos[2][k][iw]-pos[1][k][iw];//x2-x1, y2-y1 i z2-z1
                x[1][3][k]=pos[3][k][iw]-pos[1][k][iw];
                x[2][3][k]=pos[3][k][iw]-pos[2][k][iw];
                x[2][1][k]=pos[1][k][iw]-pos[2][k][iw];
                x[3][1][k]=pos[1][k][iw]-pos[3][k][iw];
                x[3][2][k]=pos[2][k][iw]-pos[3][k][iw];

                r[1][2]+=pow(x[1][2][k],2);    //r12=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2) i analogno za r23 i r13
                r[1][3]+=pow(x[1][3][k],2);
                r[2][3]+=pow(x[2][3][k],2);
                r[2][1]+=pow(x[2][1][k],2);
                r[3][1]+=pow(x[3][1][k],2);
                r[3][2]+=pow(x[3][2][k],2);
			 }
                //magnituda duljina vektora r12, r23 i r13.
				 r[1][2]=sqrt(r[1][2]);
				 r[1][3]=sqrt(r[1][3]);
				 r[2][3]=sqrt(r[2][3]);
                 r[2][1]=sqrt(r[2][1]);
				 r[3][1]=sqrt(r[3][1]);
				 r[3][2]=sqrt(r[3][2]);





		P[iw] = pow(Psi(r[1][2])*Psi(r[1][3])*Psi(r[2][3]),2);//vjerojatnost
		E[iw] = energija(r,x);//energija
	}
	SbE = 0.;
	SbE2 = 0.;
	for (ib = 1; ib <= Nb; ib++) // blokovi
	{
	    SsE = 0.;
		for (is = 1; is <= Ns; is++) // koraci
		{
			SwE = 0;
          for (iw = 1; iw <= Nw; iw++)//setaci
            {
                rnew[1][2] =0.;
                rnew[1][3] =0.;
                rnew[2][3] =0.;
				rnew[2][1] =0.;
		        rnew[3][1] =0.;
		        rnew[3][2] =0.;

                for(a=1;a<=3;a++){// petlje po parovima atomima ukupno ih je 3. r12, r13 i r23

                        for(k=1;k<=3;k++){
                            newpos[a][k]=pos[a][k][iw]+dk[k]*2.0*(ran1(&idum)-0.5);
                }}

                for(k=1;k<=3;k++){

                //međučestične komponenete.INDEX x[a][d][k] gdje je [a] = prvi atom [d]= drugi atom, [k]-komponenta x,y,z.

                xp[1][2][k]=newpos[2][k]-newpos[1][k];//x2-x1, y2-y1 i z2-z1
                xp[1][3][k]=newpos[3][k]-newpos[1][k];
                xp[2][3][k]=newpos[3][k]-newpos[2][k];
                xp[2][1][k]=newpos[1][k]-newpos[2][k];
                xp[3][1][k]=newpos[1][k]-newpos[3][k];
                xp[3][2][k]=newpos[2][k]-newpos[3][k];

                rnew[1][2]+=pow(xp[1][2][k],2);    //r12=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2) i analogno za r23 i r13
                rnew[1][3]+=pow(xp[1][3][k],2);
                rnew[2][3]+=pow(xp[2][3][k],2);
                rnew[2][1]+=pow(xp[2][1][k],2);
                rnew[3][1]+=pow(xp[3][1][k],2);
                rnew[3][2]+=pow(xp[3][2][k],2);

			 }
                //magnituda duljina vektora r12, r23 i r13.
				 rnew[1][2]=sqrt(rnew[1][2]);
				 rnew[1][3]=sqrt(rnew[1][3]);
				 rnew[2][3]=sqrt(rnew[2][3]);
                 rnew[2][1]=sqrt(rnew[2][1]);
				 rnew[3][1]=sqrt(rnew[3][1]);
				 rnew[3][2]=sqrt(rnew[3][2]);

                  Pp=pow(Psi(rnew[1][2])*Psi(rnew[1][3])*Psi(rnew[2][3]),2);//vjerojatnost u probnom pol
			      T=Pp/P[iw];//prijelaz

				// Metropolis algoritam
				if (T >= 1)
				{
				    for(a=1; a<=3;a++){
                            for (k = 1; k <= 3; k++){
					          pos[a][k][iw]=newpos[a][k];

                            }}

					accept += 1.;
					P[iw] = Pp;
					E[iw] = energija(rnew,xp);
				}
				else if (ran1(&idum) <= T)
				{
				    for(a=1; a<=3;a++){
                            for (k = 1; k <= 3; k++){
					          pos[a][k][iw]=newpos[a][k];

                            }
				    }

					accept += 1.;
					P[iw] = Pp;
					E[iw] = energija(rnew,xp);
				}

				SwE = SwE + E[iw];
			}// setaci
			if(is%30==0){
            acc_is = accept / ((ib-1)*Nw*Ns+is*Nw);
            if(acc_is > 0.5){
               for(k=1; k<=3; k++){dk[k] = dk[k]*1.05;}
            }
            if(acc_is < 0.5){
                for(k=1; k<=3; k++){dk[k] = dk[k]*0.95;}
            }}
            // akumulacija podataka nakon stabilizacije

			if (ib > NbSkip)
			{
				SsE+= SwE/Nw;
			}
		} // koraci

        if (ib > NbSkip) // akumulacija podataka nakon stabilizacije
		{
			SbE+= SsE/Ns;
			SbE2+= SsE*SsE/(Ns*Ns);

			fprintf(fout, "%7d %16.8e %16.8e\n", ib-NbSkip, SsE/Ns, SbE /(ib-NbSkip));
		}
		itmp = (int)(round(acc_is * 100.));

		 printf("%6d. blok:  %d%% prihvacenih,  Eb = %10.2e\n", ib-NbSkip, itmp, SsE / Ns);
    }// blokovi

	AE = SbE / (Nb-NbSkip);
	sigmaE = sqrt(fabs((SbE2/(Nb-NbSkip) - AE * AE)/((Nb-NbSkip)-1.)));

	E_ret[0] = AE;
	E_ret[1] = sigmaE;
	accept = accept/(Nw*Ns*Nb);
    printf("postotak prihvacenih koraka: %4.1f\n" ,accept*100.);
	printf("\n konacni max. koraci: %6.2e %6.2e %6.2e\n", dk[1], dk[2], dk[3]);
	printf("\n E = %8.5e +- %6.2e \n\n", AE, sigmaE);
    fclose(fout);
}
