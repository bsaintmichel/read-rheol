## ReadRheology 

Brice Saint-Michel, laboratoire Navier.

Ce programme vous simplifie la vie pour importer des données des rhéomètres Anton Paar RheoCompass, Malvern rSpace et TA Trios, du format .csv ou .txt vers Python, puis de les tracer simplement à l'aide du paquet [Matplotlib](https://matplotlib.org/), facilement installable avec `pip` (et probablement déjà sur votre distribution Python). 

Trois fichiers d'exemple (un d'Anton Paar, un de Malvern, un de TA) sont fournis à titre d'exemple.

------------------------------

### Programme de base

Un exemple d'utilisation du programme se trouve dans le fichier [ReadRheology.ipynb](./ReadRheology.ipynb), section "Utilisation de base"

Le fichier appelle des fonctions du fichier Python `rheol_functions.py` pour fonctionner, notamment : 
 
* le code principal qui lit les fichiers `.csv` : `read_rheology()`. Le programme détecte les fichiers Malvern ou Anton Paar, découpe proprement les étapes (step / interval), gère (normalement) les données de LAOS et renvoie un tableau (Pandas Dataframe). Il renomme les colonnes courantes (contrainte, déformation, etc.) en des identifiants simples (`stress`, `strain`, ...), et conservera les noms plus compliqués des colonnes qu'il ne reconnaît pas automatiquement.
* une fonction de 'confort' :  `slice` permettant de renvoyer un sous-tableau contenant les étapes qui vous intéressent, choisies par numéro d'étape ou par type (voir ci-dessous).
* une autre fonction de confort :  `assign_steps`, permettant de définir un 'type' pour chaque étape, qui pourra être reconnu plus tard par les fonctions de tracé et les fonctions de confort.
* des fonctions de tracé :  `plot_flowcurve()`, `plot_fsweep()`, `plot_asweep()`, `plot_tsweep()` qui permettent de tracer correctement vos courbes d'écoulement, balayages en amplitude et en fréquence, etc. Il est bien sûr possible de superposer ces tracés si vous leur donnez à manger un tableau contenant plusieurs étapes.

--------

### Comment exporter ses données au format .csv ?

1. Anton Paar RheoCompass

Débrouillez-vous pour afficher un tableau avec les quantités physiques qui vous intéressent (vous pouvez le faire en définissant vous-même un tableau modèle dans votre projet, ou sinon vous pouvez le faire à la main à chaque fois). 

Dans la vue "tableau", cochez toutes les séries de données que vous voulez afficher. Ensuite, allez dans "Accueil" et cliquez sur "Exporter l'objet". Il est possible que cela ratouille un peu, et qu'on vous propose d'enregistrer au format `.RhPts`. Dans ce cas, jouez à cliquer sur des objets différents (une étape, i.e. les vignettes d'oscilloscope plutôt qu'un test, i.e. la boîte en carton) sans rien cocher/décocher (juste changer de sélection) et normalement vous pourrez finir par exporter en `.csv` . 

2. Malvern rSpace

Si vous êtes à Navier, vous devrez créer une "table" qui contient les résultats, et utiliser le modèle de tableau "NAVIER_DEFAULT" qui devrait normalement contenir toutes les valeurs qui comptent. Je crois que le modèle rajoute automatiquement les données de LAOS quand vous "glissez" les tests du menu de gauche vers le tableau pour les ajouter. 

Cliquez ensuite (clic droit) sur le coin en haut à gauche du tableau, et faites 'export --> to .csv' et le tour est joué

3. Trios

Sachez que vous pouvez cocher pour afficher des "variables" supplémentaires lorsque vous sélectionnez vos variables dans la vue tableau ou graphe : c'est une petite case à cocher en bas à droite de la fenêtre de dialogue.

Quand vous avez sélectionné vos séries de données, faites un clic-droit sur le fichier source lui-même, ou cliquez sur le logo TA en haut à gauche du menu, et faites "Exporter" au format "texte" (sans délimitation par des virgules). Cochez la case correspondant aux paramètres de l'expérience (qui rajoute des lignes correspondant aux détails de la géométrie parfois nécessaires). Et récupérez votre fichier .txt.

-----------------------

### Utilisation plus avancée des fonctions de LAOS (Anton Paar / Malvern): 

Un exemple d'utilisation du programme se trouve dans le fichier [ReadRheology.ipynb](./ReadRheology.ipynb), section "LAOS (Anton Paar)" ou "LAOS (Malvern)" en fonction du rhéomètre utilisé. Cette fois-ci, on va utiliser en plus des fonctions précédente deux petits outils utiles : 

* `proj_fourier()` :  qui projette un signal sur une série de Fourier. Il est un peu bête, donc il faut lui donner à manger un temps normalisé par la fréquence angulaire (0 à 2 $\pi$ pour une période) et c'est mieux de lui donner un nombre entier et pair de périodes. Sinon vos coefficients de la projection vont déconner. Le programme renvoie un `dict` contenant à la fois la partie réelle et imaginaire des harmoniques $k$, (respectivement `sin` et `cos`), mais également leur module (`amp`) et phase (`phs`). 
* `build_fourier()` :  qui reconstruit le signal à partir du dictionnaire des harmoniques $k$ et une base de temps, pour notamment vérifier la qualité de votre projection.