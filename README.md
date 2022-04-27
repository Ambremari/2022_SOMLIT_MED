# Caractériser la structure et la dynamique des assemblages phytoplanctoniques du milieu côtier de la série temporelle SOMLIT
## Approches statistiques appliquées aux données fonctionnelles de la période 2012-2021

Vous trouverez ici les données et les codes permettant de reproduire les analyses réalisées durant mon stage de fin de Master 2 de janvier à juin 2022. 
L'analyse est réalisée sur les données hydrologiques et les données de cytométrie du pico et nano-plancton pour les trois stations de Méditerrannée 
faisant partie du réseau SOMLIT (Banyuls, Marseille et Villefranche). 

### Code

Preprocessing.R permet d'extraire et de mettre en forme les données nécessaires aux analyses statistiques

Regression.R réalise les régressions quadratiques locales avec noyau gaussien pour les séries temporelles des données brutes, avec l'intervalle de confiance à 95% obtenus à partir d'une approximation de Taylor d'ordre 2, et construit les figures associées. 

Decomposition.R permet de réaliser une décomposition MSTL des régressions, de supprimer la corrélation linéaire des résidus induite par suréchantillonnage et de visualiser le périodogramme de la série, la série décomposée, la corrélation linéaire des résidus, les résidus décorrélés, l'ACF des résidus et l'histogramme des résidus.  

### Données

Les données sont issues du site du Service d’Observation en Milieu Littoral (https://www.somlit.fr/)

HYDRO.csv : données hydrologiques SOMLIT de 1994 à 2021 pour les trois stations de Méditerranée : Banyuls, Marseille et Villefranche. 

PICONANO.csv : données du pico et nano-plancton SOMLIT de 2012 à 2021 pour les trois stations de Méditerranée : Banyuls, Marseille et Villefranche. 

CTD_MARSEILLE.csv : données CTD SOMLIT de 1994 à 2021 pour Marseille. 

CTD_BANYULS.csv : données CTD SOMLIT de 1994 à 2021 pour Banyuls.

CTD_VILLEFRANCHE.csv : données CTD SOMLIT de 1994 à 2021 pour Villefranche. 

PICONANO_AB.csv : données d'abondance du pico et nano-plancton SOMLIT de 2012 à 2021 pour les trois stations de Méditerranée : Banyuls, Marseille et Villefranche.

PICONANO_DIFF.csv : données de diffusion lumineuse du pico et nano-plancton SOMLIT de 2012 à 2021 pour les trois stations de Méditerranée : Banyuls, Marseille et Villefranche.

### Contact 

email : m.couteyen_carpaye@ecologie.re
