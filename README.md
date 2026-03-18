Ce dépôt git alimente le projet EMI d'Arthur DEROISSART, Samuel BESLON et Paul BUET.

Il a pour but de modéliser la sédimentation d'abord dans un modèle 1D à travers différents schémas. En vue de pouvoir les comparer aux schémas 
actuellement utilisés dans les modèles opérationnels 3D (type mésoNH ou AROME).

# Projet EMI 2026 Phy-nh : "Vers une revisite de la sédimentation dans Arome"

Le projet vise à développer différents alrorythmes implémentant le schéma de sédimentation Box Lagrangien en utilisant différentes méthodes d'avance temporelle
(Step by step ou step forward) et en faisant ou non varier la taille de la boite utilisée pour ce schéma. Les resultats obtenus pourront être comparés à des
modèles utilisant le schéma eulérien utilisé en opérationnel.

## Installations et configurations nécessaires

<ul>
<li> pour se mettre sur l'environnement : python -m venv mon_env ; source mon_env/bin/activate ; pip install -r requirements.txt
<li> pour configurer le branchement phyex : cloner le dépôt git PHYEX en amont du dépôt git de ce projet puis executer activation branchement.sh dans le dépôt de ce projet.
<li> pour lancer les scripts : python3 Scripts/entree.py -h
<li> le module Argparse vous renverra les informations pour paramétrer votre run.

### Architecture du projet

L'appel des scripts se fait via le script entree.py lui-même utilisant le module Argparse pour recevoir les informations et les communiquer à distribution.py
chargé de l'appel des modèles en fonction des informations reçues de l'utilisateur.

Une fois les informations fournies à distribution, ce script appelle le script du modèle correspondant. Parmi ces scripts on retrouve :

<ul><li> <br> box_lagrangien.py </br> (modèle par défaut à largeur de boite fixe et utilisant la méthode d'avance temporelle <br> Step_By_Step </br>)
<li> <br> box_lagrangien_vectorized.py </br> (même modèle que ci-dessus mais incluant de légères modifications concernant la vectorisation des porcédés)
<li> <br> box_lagrangien_sf.py </br> (modèle à largeur de boite fixe utilisant la méthode d'avance temporelle <br>Step_Forward</br>)
<li> <br> box_lagrangien_sf_vectorized.py </br> (même modèle que ci-dessus mais intégrant de légères modifications concernant la vectorisation des porcédés)
<li> <br> phyex.py </br> (code d'interfaçage permettant de lancer les modèles <br>EULE, EULE2 et STAT </br>en utilisant les mêmes formats de données et les mêmes conditions initiales)
<li> <br> box_lagrangien_def.py </br> (modèle à largeur de boite variable utilisant la méthode d'avance temporelle <br>Step_By_Step</br>)
