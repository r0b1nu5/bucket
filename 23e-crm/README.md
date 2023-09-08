# CRM - Colloque 2023 - Champéry
## Graphes et réseaux électriques

*En cas de questions : robin.delabays@hevs.ch*

Vous trouverez ici quelques codes en accompagnement de la présentation (le pdf de la présentation sera mis à disposition ici).

### Prérequis
Les codes sont écrits en Julia, un langage similaire an Matlab, mais en accès libre. L'installation se fait directement depuis [le site de Julia](https://julialang.org). 

En plus des librairies de base, nous allons utiliser les librairies spécifiques suivantes :

- DelimitedFiles
- Graphs
- LinearAlgebra
- PyPlot
- Statistics

Vous pouvez copier le code suivant dans le terminal Julia pour installer les librairies nécessaires :

`using Pkg; Pkg.add("DelimitedFiles"); Pkg.add("Graphs"); Pkg.add("LinearAlgebra"); Pkg.add("PyPlot"); Pkg.add("Statistics");`

Il faut ensuite télécharger l'ensemble des fichier de ce répertoire, soit en clonant le lien GitHub, soit en téléchargeant et décompressant le fichier compressé `23e-crm.zip`. 

### Lancer un script
Il y a deux scripts à disposition :

- `crm-ex-nodes.jl`
- `crm-ex-lines.jl`

Pour lancer un script, il faut commencer par lancer Julia, ce qui ouvre un terminal. 
De là, il faut se rendre dans le répertoire où les codes et données ont été téléchargées :

`cd("/Adresse/Du/Répertoire/")`

On lance ensuite le script par l'instruction :

`include("crm-ex-nodes.jl")` ou `include("crm-ex-lines.jl")`

La suite des instructions s'affichera dans le terminal.


