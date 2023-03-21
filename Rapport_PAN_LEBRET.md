## Parallélisation du code

### Parallélisation en mémoire partagée

On utilisera la simulation du model avec 3 vortexs pour observer l'accélération du a OPEN MP.

Sans OPEN MP ,on tourne aux alentours de 15 FPS.
Avec OPEN MP ,on parallélise la première boucle dans solve_RK4_movable_vortices et on voit une augmentation de 15 fps.
En parallélisant la deuxième boucle on voit une nouvelle augmentation d'une 30 aine de FPS et on passe à 60

Cependant avec la partie suivante, la parallélisation de la seconde boucle à l'air de réduire considérablement la parallélisation de la seconde partie.

### Parallélisation en mémoire distribuée et partagée des calculs

Pour la suite, la parallélisation avec plus de 2 process est dans la branche grid separation

L'idée est simple, le processus 0 s'occupe de l'affichage, le processus 1 s'occupe de communiquer avec le processus 0 et avec les autres processus de calcul. //
On crée un communicateur contenant tout les processus faisant du calcul pour pouvoir faire des broadcast sans encombrer le processus 0
Le programme à l'air de fonctionner correctement pour n = 4 ou 8
