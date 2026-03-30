Ce dépôt git alimente le projet EMI d'Arthur DEROISSART, Samuel BESLON et Paul BUET.

Il a pour but de modéliser la sédimentation d'abord dans un modèle 1D à travers différents schémas. En vue de pouvoir les comparer aux schémas 
actuellement utilisés dans les modèles opérationnels 3D (type mésoNH ou AROME).

# Projet EMI 2026 Phy-nh : "Vers une revisite de la sédimentation dans Arome"

Le projet vise à développer différents alg  orythmes implémentant le schéma de sédimentation Box Lagrangien en utilisant différentes méthodes d'avance temporelle
(Step by step ou step forward) et en faisant ou non varier la taille de la boite utilisée pour ce schéma. Les resultats obtenus pourront être comparés à des
modèles utilisant le schéma eulérien utilisé en opérationnel.

## Installations et configurations nécessaires

<ul>
<li> pour se mettre sur l'environnement : python -m venv mon_env ; source mon_env/bin/activate ; pip install -r requirements.txt
<li> pour configurer le branchement phyex : cloner le dépôt git PHYEX en amont du dépôt git de ce projet puis executer activation branchement.sh dans le dépôt de ce projet.
<li> pour lancer les scripts : python3 Scripts/entree.py -h
<li> le module Argparse vous renverra les informations pour paramétrer votre run.
</ul>

## Architecture du projet
### Fonctionnement général des scripts

Les scripts python utilisés sont tous contenus dans "scripts".

L'appel des scripts se fait via le script entree.py lui-même utilisant le module Argparse pour recevoir les informations et les communiquer à distribution.py
chargé de l'appel des modèles en fonction des informations reçues de l'utilisateur.

Une fois les informations fournies à distribution, ce script appelle le script du modèle correspondant. Parmi ces scripts on retrouve :

<ul>
<li>
box_lagrangien.py (modèle par défaut à largeur de boite fixe et utilisant la méthode d'avance temporelle Step_By_Step)
</li>
<li>
box_lagrangien_vectorized.py (même modèle que ci-dessus mais incluant de légères modifications concernant la vectorisation des porcédés)
</li>
<li>
box_lagrangien_sf.py (modèle à largeur de boite fixe utilisant la méthode d'avance temporelle Step_Forward)
</li>
<li>
box_lagrangien_sf_vectorized.py (même modèle que ci-dessus mais intégrant de légères modifications concernant la vectorisation des porcédés)
</li>
<li>
phyex.py (code d'interfaçage permettant de lancer les modèles EULE, EULE2 et STAT en utilisant les mêmes formats de données et les mêmes conditions initiales)
</li>
<li>
box_lagrangien_def.py (modèle à largeur de boite variable utilisant la méthode d'avance temporelle Step_By_Step)
</li>
</ul>


Les modèles appellent tous condi_init.py pour initialiser les prodils de contenu et/ou de masse (selon l'usage d'un schéma bin ou d'un schéma bulk).

Les resultats des modèles sont ensuite récupérés dans distibution.py.

Enfin celui-ci fait appel à affichage.py pour afficher les resultats.

### Description du programme entrée

C'est le premier fichier appelé par le modèle, le programme entrée.py permet à l'utilisateur de spécifier des conditions spécifiques utilisé ensuite par le modèle comme par exemple la durée de la simulation ou encore le nombre de bin à utiliser dans le cas d'un modèle Box_Lagrangien.

Ce fichier passe en revu chacun des paramètres spécifiques entrée par l'utilisateur pour verifier qu'il ait été rentré sous la bonne forme (nombre entier si le paramètre doit être un nombre entier, "Yes" or "No" s'il s'agit d'un choix...). Une fois que chacun des paramètre est vérifié, le fichier lance le programme distribution.

### Description du programme distribution

Ce prgogramme permet d' éxecuter le bon fichier en fonction du modèle choisit par l'utilisateur. Par exemple, si l'utilisateur choisit d'utiliser un modèle Box Lagrangien, déformable en Step By Step, le programme distribution va executer le fichier box_lagrangien_def.py.

Chaque modèle renvoie un résultat, stocké par la suite dans la variable result. Il se compose généralement de 3 listes : une première qui correspond à l'évolution dans le temps de la concentration dans chacune des mailles (liste de liste), une seconde qui correspond à la masse de précipitation au sol (liste simple) et une troisième qui correspond à l'évolution dans le temps de la masse dans chacune des mailles (liste de liste).

Une fois les modèles executés, le programme distribution lance le calcul de la courbe théorique de cumul qu'on est censé obtenir. Puis appel le fichier affichage qui s'occupe d'afficher et d'enregistrer les graphes.

### Description du programme d'affichage

Le programme commence par une courte initialisation qui permet de changer le chemin d'enregistrement en fonction du modèle appelé par l'utilisateur. Si il s'agit d'un modèle déformable, il n'y a pas besoin de renseigner le nombre de bin. De mème s'il s'agit d'un modèle Statistique ou Eulerien (1 ou 2 moment), il n'y a pas besoin de dire si le modèle est deformable ou de spécifier le type d'avance temporelle.

La variable params permet de stocker l'ensemble des paramètres modèles qu'il est utile de mentionner sur les figures. L'objectif étant de savoir à quoi correspond un graphe quand on l'analyse.

Le programme execute ensuite l'affichage des graphes. L'ajout de plt.close() à la fin de chacune des deux fonctions (une fonction pour la sédimentation et une fonction pour les précipitations au sol), permet de lancer plusieurs runs à la suite. Sans cette ligne, le modèle s'arrête tant que les graphes sont ouverts et poursuit une fois que les graphes sont fermés manuellement. Les figures sont enregistrées dans les fichiers correspondant aux paramètres entrée initalement par l'utilisateur.

### Stockage des figures dans "fig"

Le dossier fig contient les figures produites par les scripts. Il existe cependant une option à spécifier au moment de l'appel permettant de les produire dans un dossier "fig" situé ailleurs. 
Les fichier seront enregistrés dans un chemin sous la forme nom_du_modèle/paramètre_suivi.png.
Attention, il est nécessaire pour leur enregistrement que les répértoires nom_du_modèle soient créés avant l'execution des modèles.

### Documents annexes

