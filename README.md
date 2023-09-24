# Convertisseur AIXM vers KML

Ce script Python est une adaptation de [aixm2Parser](https://github.com/cquest/aixmParser) il est utilisé pour convertir les données AIXM 4.5(Aeronautical Information Exchange Model) en format KML (Keyhole Markup Language) qui est utilisé pour la visualisation géospatiale.

## Dépendances

Ce script utilise plusieurs bibliothèques Python, notamment:
- `simplekml` pour générer des fichiers KML.
- `shapely` pour traiter et manipuler les données géospatiales.
- `BeautifulSoup` de `bs4` pour parser les fichiers XML.
- `pyproj` pour la projection des coordonnées.

Assurez-vous d'installer ces dépendances avant d'exécuter le script.

## Comment utiliser

### Arguments de la ligne de commande

Ce script supporte plusieurs arguments de la ligne de commande pour sélectionner les types de données à exporter :
- `--ahp` : Exporter les aérodromes/héliports.
- `--obs` : Exporter les obstacles.
- `--rcp` : Exporter les centres de piste.
- `--gbr` : Exporter les frontières géographiques.
- `--abd` : Exporter les limites de l'espace aérien.
- `--uni` : Exporter les tours de contrôle.
- `--gsd` : Exporter les emplacements des portes.

### Exemple d'exécution

Pour convertir un fichier AIXM en KML, utilisez la commande suivante :
```
./nom_du_script.py --ahp --obs chemin_vers_fichier_AIXM.aixm
```

Ceci convertira les aérodromes/héliports et les obstacles du fichier AIXM spécifié.

## Fonctionnalités

Ce script offre plusieurs fonctionnalités, telles que :
- La conversion de différents types de données géospatiales.
- La génération de couleurs aléatoires pour la visualisation.
- La prise en charge des polygones, des lignes et des points.

NB: Pour l'heure le script est en V0.1 