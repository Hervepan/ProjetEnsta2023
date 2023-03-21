## Parallélisation du code

### Parallélisation en mémoire partagée

Sans OPEN MP en se basant sur la modélisation des 3 vortex, on tourne aux alentours de 15 FPS.
Avec OPEN MP on parallélise la boucle dans solve_RK4_movable_vortices et on voit une augmentation de 15 fps. 
Cependant avec la partie suivante, OPEN MP a l'air de réduire les performances surement du aux grand nombres de communication 
entre les process

### Parallélisation en mémoire distribuée et partagée des calculs

L'idée est simple, le processus 0 s'occupe de l'affiche, le processus 1 est le maitre des processus de calculs. //
Il distribue la charge de travail aux autres processus et ensuite récupère les nouvelles données.
