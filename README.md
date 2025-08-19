## 2. Varijacijski Monte Carlo

Razmotrimo homogeni sustav identičnih čestica koje međudjeluju s uparenim isključivo radijalnim međučestičnim potencijalom $V(r_{ij})$. Hamiltonijan je dan:

$$
\mathcal{H} = -\frac{\hbar^2}{2m}\sum_{i=1}^{N} \nabla_i^2 + \sum_{i<j} V(r_{ij}) = \mathcal{T} + \mathcal{V}. \tag{2.1}
$$

Varijacijska metoda polazi od probne valne funkcije $\psi(\vec{r},\vec{a})$ koja ovisi o varijacijskim parametrima i položajima atoma u sustavu. Cilj je naći parametre koji minimiziraju funkcional energije:

$$
E = \frac{\langle \psi | H | \psi \rangle}{\langle \psi | \psi \rangle} \geq E_0. \tag{2.2}
$$

Varijacijski princip iskazuje da je za bilo koju valnu funkciju očekivana vrijednost od Hamiltonijana viša od osnovnog stanja egzaktne energije $E_0$.  
Izračun energije za danu probnu valnu funkciju je nedeterministički zadatak jer se treba izračunati višedimenzionalni integral:

$$
E = \frac{\int d^3r_1 \ldots d^3r_N \, \psi^*(\vec{r},\vec{a}) H \psi(\vec{r},\vec{a})}{\int d^3r_1 \ldots d^3r_N |\psi(\vec{r},\vec{a})|^2}. \tag{2.3}
$$

Kako je poznato, višedimenzionalna integracija je relativno lagana za Monte Carlo. Varijacijskim Monte Carlom je moguće izračunati egzaktnu energiju do na statistički šum.  
Definirajmo funkciju gustoće vjerojatnosti:

$$
f(R) = \frac{|\psi(R)|^2}{\int dR |\psi(R)|^2}, \tag{2.4}
$$

a lokalna energija je

$$
E_L(R) = \frac{1}{\psi(R)} H \psi(R), \tag{2.5}
$$

očekivana vrijednost Hamiltonijana je dana

$$
\langle H \rangle_\psi = \int dR \, E_L(R) f(R). \tag{2.6}
$$

Energija se dobiva kao srednja vrijednost od $E_L(R)$ pri čemu se položaji šetača/konfiguracija atoma $R = r_1 \ldots r_n$ generiraju sukladno iz distribucije valne funkcije $f(R)$:

$$
\langle H \rangle_\psi = \frac{1}{n} \sum_{i=1}^{n} E_L(R_i). \tag{2.7}
$$

gdje je $n$ broj šetača. Varijanca se određuje po formuli gdje je $n_b$ broj blokova:

$$
\sigma^2 = \frac{1}{n_b - 1}\left[ \frac{1}{n_b} \sum_{i=1}^{n_b} E_L^2(R_i) - \left( \frac{1}{n_b} \sum_{i=1}^{n_b} E_L(R_i) \right)^2 \right]. \tag{2.8}
$$

Vrijednost funkcionala (2.7) računamo Monte Carlo integracijom, pri čemu se koristi **Metropolisov algoritam** koji prihvaća one položaje/pomake šetača tako da njihova distribucija odgovara valnoj funkciji vjerojatnosti $f(R)$.  

Koraci algoritma su sljedeći:

1. Inicijalizirati šetače na početne položaje $R_i^0 \ (i=1,\ldots,n)$  
2. Predložiti pomake $R_i^f = R_i^0 + (2 \cdot (\text{ran(idum)} - 0.5))$  
3. Izračunamo prijelaz:  

   $$
   T(R_i^0 \to R_i^f) = \frac{|\psi(R_i^f)|^2}{|\psi(R_i^0)|^2}. \tag{2.9}
   $$

4. Ako je $T > 1$ prihvaćamo pomak i $R_i^0 = R_i^f$  
5. Inače ako je $T < 1$, generiramo slučajan broj iz intervala $[0,1]$ i ako je slučajan broj manji od $T$ onda $R_i^0 = R_i^f$  
6. Na taj način uzorkujemo položaje čija distribucija je jednaka vjerojatnosti valne funkcije $|\psi|^2$ i onda računamo energiju kao prosjek po šetačima, koracima i blokovima  

---

Lokalnu energiju možemo zapisati kao zbroj kinetičkog i potencijalnog dijela:

$$
E_L(R) = \frac{(\mathcal{T} + \mathcal{V}) \psi(\vec{r},\vec{a})}{\psi(\vec{r},\vec{a})} = E_L^T(R) + E_L^V(R). \tag{2.10}
$$

gdje su $\mathcal{T}$ i $\mathcal{V}$ operatori iz (2.1).  
Valnu funkciju zapisujemo kao produkt dvočestičnih korelacijskih funkcija:

$$
\psi(R) = \prod_{i<j}^{N_{atom}} f_{ij}(r_{ij}). \tag{2.11}
$$

Doprinos kinetičkog dijela se može zapisati na optimiziran način:

$$
E_L^T(R) = \sum_{i=1}^N D_i \left[ \sum_{j \neq i}^N f_{ij}^{dr}(r_{ij}) \left( \sum_{k \neq i}^N f_{ik}^{dr}(r_{ik}) + \frac{2}{r_{ij}} \right) + \sum_{j \neq i}^N f_{ij}^{ddr}(r_{ij}) \right], \tag{2.12}
$$

gdje su $D_i = \frac{\hbar^2}{2m}$, a funkcije su:

$$
f_{ij}^{dr}(r_{ij}) = \frac{df(r)}{dr} \cdot \frac{1}{f(r) \cdot r}, \tag{2.13}
$$

i

$$
f_{ij}^{ddr}(r_{ij}) = \frac{d f_{ij}^{dr}(r_{ij})}{dr} \cdot r + 3 f_{ij}^{dr}(r_{ij}). \tag{2.14}
$$