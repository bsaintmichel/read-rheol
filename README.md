### ReadRheology 

Brice Saint-Michel, laboratoire Navier.

Ce programme vous simplifie la vie pour importer des données des rhéomètres Anton Paar [RheoCompass] et Malvern [rSpace], du format .csv ou .txt vers Python, puis de les tracer simplement à l'aide du paquet [Bokeh](https://bokeh.org/), facilement installable avec `pip`. 


------------------------------

#### Programme de base

Un exemple d'utilisation du programme se trouve dans le fichier [ReadRheology.ipynb](./ReadRheology.ipynb), section "Utilisation de base"

Le fichier appelle des fonctions du fichier Python `rheol_functions.py` pour fonctionner, notamment : 
 
* le code principal qui lit les fichiers `.csv` : `read_rheology()`. Le programme détecte les fichiers Malvern ou Anton Paar, découpe proprement les étapes (step / interval), gère (normalement) les données de LAOS et renvoie un tableau (Pandas Dataframe). Il renomme les colonnes courantes (contrainte, déformation, etc.) en des identifiants simples (`stress`, `strain`, ...)
* une fonction de 'confort' :  `slice` permettant de renvoyer un sous-tableau contenant les étapes qui vous intéressent.
* une autre fonction de confort :  `assign_steps`, permettant de définir un 'type' pour chaque étape, qui pourra être reconnu plus tard par les fonctions de tracé.
* des fonctions de tracé :  `plot_flowcurve()`, `plot_fsweep()`, `plot_asweep()`, `plot_tsweep()` qui permettent de tracer correctement vos courbes d'écoulement, balayages en amplitude et en fréquence, etc. Il est bien sûr possible de superposer ces tracés si vous leur donnez à manger un tableau contenant plusieurs étapes.

--------

#### Utilisation plus avancée des fonctions de LAOS : 

Un exemple d'utilisation du programme se trouve dans le fichier [ReadRheology.ipynb](./ReadRheology.ipynb), section "Utilisation avancée en LAOS". Cette fois-ci, on va utiliser en plus des fonctions précédente deux petits outils utiles : 

* `proj_fourier()` :  qui projette un signal sur une série de Fourier. Il est un peu bête, donc il faut lui donner à manger un temps normalisé par la fréquence angulaire (0 à 2 $\pi$ pour une période) et c'est mieux de lui donner un nombre entier et pair de périodes. Sinon vos coefficients de la projection vont déconner. Le programme renvoie un `dict` contenant à la fois la partie réelle et imaginaire des harmoniques $k$, (respectivement `sin` et `cos`), mais également leur module (`amp`) et phase (`phs`). 
* `build_fourier()` :  qui reconstruit le signal à partir du dictionnaire des harmoniques $k$ et une base de temps, pour notamment vérifier la qualité de votre projection.