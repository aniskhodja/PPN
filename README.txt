-------------------------------EXCECUTION DU PROGRAMME-----------------------------------

>>	make

	un excécutable (run) est généré 

*****************************PRODUIT SCALAIRE********************************************

>>	mpirun -np p ./run m

	remplacez les arguments p et m par des chiffre 
	p: nombre de processus
	m: taille des deux vecteurs



*****************************PRODUIT MATRICE VECTEUR*************************************

>>	mpirun -np p ./run m n k

	remplacez les arguments p , m , n , k par des chiffre 
	p: nombre de processus
	m: nombre de ligne matrice
	n: nombre de colonne matrice
	k: taille du vecteur, doit etre egale à n
	
	

*****************************PRODUIT MATRICE MATRICE*************************************

>>	mpirun -np p ./run m n k

	remplacez les arguments p , m , n , l , k par des chiffre 
	p: nombre de processus
	m: nombre de ligne matrice 1
	n: nombre de colonne matrice 1
	l: nombre de ligne matrice 2 , doit etre egale à n
	k: nombre de colonne matrice 2
