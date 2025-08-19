# Monte-Carlo-simulacije-helija3
Energije i valne fje molekule trimera He_3 .
\section{Varijacijski Monte Carlo}

Razmotrimo homogeni sustav identičnih čestica koje međudjeluju s uparenim isključivo radijalnim međučestičnim potencijalom $V(r_{ij})$. Hamiltonijan je dan \cite{2}:

\begin{equation}
\mathcal{H} = -\frac{\hbar^2}{2m}\sum_{i=1}^{N} \nabla_i^2 + \sum_{i<j} V(r_{ij}) = \mathcal{T} + \mathcal{V}.
\tag{2.1}
\end{equation}

Varijacijska metoda polazi od probne valne funkcije $\psi(\vec{r},\vec{a})$ koja ovisi o varijacijskim parametrima i položajima atoma u sustavu. Cilj je naći parametre koji minimiziraju funkcional energije

\begin{equation}
E = \frac{\langle \psi | H | \psi \rangle}{\langle \psi | \psi \rangle} \geq E_0.
\tag{2.2}
\end{equation}

Varijacijski princip iskazuje da je za bilo koju valnu funkciju očekivana vrijednost od Hamiltonijana je viša od osnovnog stanja egzaktne energije $E_0$. Izračun energije za danu probnu valnu funkciju je nedeterministički zadatak jer se treba izračunati višedimenzionalni integral

\begin{equation}
E = \frac{\int d^3r_1 \ldots d^3r_N \psi^*(\vec{r},\vec{a}) H \psi(\vec{r},\vec{a})}{\int d^3r_1 \ldots d^3r_N |\psi(\vec{r},\vec{a})|^2}.
\tag{2.3}
\end{equation}

U ovom je kontekstu teorijski problem Monte Carlo metoda vrlo korisna. Kao što je poznato, višedimenzionalna integracija je relativno lagana za Monte Carlo. Varijacijskim Monte Carlom je moguće izračunati egzaktnu energiju do na statistički šum. Definirajmo funkciju gustoće vjerojatnosti

\begin{equation}
f(R) = \frac{|\psi(R)|^2}{\int dR |\psi(R)|^2},
\tag{2.4}
\end{equation}

a lokalna energija je

\begin{equation}
E_L(R) = \frac{1}{\psi(R)} H \psi(R),
\tag{2.5}
\end{equation}

očekivana vrijednost hamiltonijana je dana

\begin{equation}
\langle H \rangle_\psi = \int dR \, E_L(R) f(R),
\tag{2.6}
\end{equation}

energija se dobiva kao srednja vrijednost od $E_L(R)$ pri čemu se položaji šetača/konfiguracija atoma 
$R = r_1 \ldots r_n$ generiraju sukladno iz distribucije valne funkcije $f(R)$

\begin{equation}
\langle H \rangle_\psi = \frac{1}{n} \sum_{i=1}^{n} E_L(R_i),
\tag{2.7}
\end{equation}

gdje je $n$ broj šetača. Varijanca se određuje po formuli gdje je $n_b$ broj blokova

\begin{equation}
\sigma^2 = \frac{1}{n_b - 1}\left[ \frac{1}{n_b} \sum_{i=1}^{n_b} E_L^2(R_i) - \left( \frac{1}{n_b} \sum_{i=1}^{n_b} E_L(R_i) \right)^2 \right].
\tag{2.8}
\end{equation}

Vrijednost funkcionala (2.7) računamo Monte Carlo integracijom, pri čemu se koristi \textbf{Metropolisov} algoritam koji prihvaća one položaje/pomake šetača tako da njihova distribucija duljina odgovara valnoj funkciji ili vjerojatnosti $f(R)$. Koraci algoritma su sljedeći:

\begin{enumerate}
\item Inicijalizirati šetače na početne položaje $R_i^0 \ (i=1,\ldots,n)$,
\item Predložiti pomake $R_i^f = R_i^0 + (2 \cdot (\text{ran(idum)} - 0.5))$,
\item Izračunamo prijelaz:
\begin{equation}
T(R_i^0 \to R_i^f) = \frac{|\psi(R_i^f)|^2}{|\psi(R_i^0)|^2}.
\tag{2.9}
\end{equation}
\item Ako je $T > 1$ prihvaćamo pomak i $R_i^0 = R_i^f$,
\item Inače ako je $T < 1$, generiramo slučajan broj iz intervala $\text{ran}() \in [0,1]$ i ako je slučajan broj manji od $T$ onda $R_i^0 = R_i^f$,
\item Na taj način uzorkujemo položaje čija distribucija je jednaka vjerojatnosti valne funkcije $|\psi|^2$ i onda računamo energiju kao prosjek po šetačima, koracima i blokovima.
\end{enumerate}

Lokalnu energiju možemo zapisati kao zbroj kinetičkog i potencijalnog dijela

\begin{equation}
E_L(R) = \frac{(\mathcal{T} + \mathcal{V}) \psi(\vec{r},\vec{a})}{\psi(\vec{r},\vec{a})} = E_L^T(R) + E_L^V(R),
\tag{2.10}
\end{equation}

gdje su $\mathcal{T}$ i $\mathcal{V}$ operatori iz (2.1).  
Valnu funkciju zapisujemo kao produkt dvočestičnih korelacijskih funkcija

\begin{equation}
\psi(R) = \prod_{i<j}^{N_{atom}} f_{ij}(r_{ij}).
\tag{2.11}
\end{equation}

Doprinos kinetičkog dijela se može zapisati na optimiziran način čiji izvod je raspisan u literaturi, a može se numerički izraziti kao

\begin{equation}
E_L^T(R) = \sum_{i=1}^N D_i \left[ \sum_{j \neq i}^N f_{ij}^{dr}(r_{ij}) \left( \sum_{k \neq i}^N f_{ik}^{dr}(r_{ik}) + \frac{2}{r_{ij}} \right) + \sum_{j \neq i}^N f_{ij}^{ddr}(r_{ij}) \right],
\tag{2.12}
\end{equation}

gdje su $D_i = \frac{\hbar^2}{2m}$, a funkcije su

\begin{equation}
f_{ij}^{dr}(r_{ij}) = \frac{df(r)}{dr} \cdot \frac{1}{f(r) \cdot r},
\tag{2.13}
\end{equation}

i

\begin{equation}
f_{ij}^{ddr}(r_{ij}) = \frac{d f_{ij}^{dr}(r_{ij})}{dr} \cdot r + 3 f_{ij}^{dr}(r_{ij}).
\tag{2.14}
\end{equation}

