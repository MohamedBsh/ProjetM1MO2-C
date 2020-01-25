/**************************************************************/
/*  Projet de Programmation Partie I                          */
/*  R�alis� par : Mohamed-Amine Bousahih                      */
/*  Num�ro etd  : 21500267                                    */    
/*  Universit� Paris Diderot                                  */
/*  Master 1 Math�matiques et Applications                    */
/*  Subject : Identification of Common Molecular Subsequences */
/**************************************************************/

/**************************** Librairies ****************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**************************** Structures ****************************/

typedef struct // structure qui contient une valeur et son emplacement dans la matrice.
{
    int idx;
    int idy;
    double value;
} cell;

typedef struct  // structure qui contient l'�l�ment maximum de la matrice ainsi que son nombre d'occurences.
{
	double val_max;
	int count;
} cell_max;

struct path  // structure qui contient une structure cell et un pointeur sur la structure suivante. Element de notre liste chain�e.
{
    cell e;
    struct path* next;
};
typedef struct path path;

typedef struct  // structure controlant notre liste chain�e.
{
    path* tete;
} cell_liste;


/***********************************        Prototypes            **************************************/


/********************************   Pour la matrice de notation   ***************************************/

double** allocation(int l1, int l2); // allocation dynamique de notre matrice de notation. 
double similarite (char s1, char s2); // renvoie un coefficient de similarit� entre deux caract�res. 
double maxQuatrecoeff (double a, double b, double c, double d); // renvoie le coefficient maximal entre 4 valeurs (selon la m�thode de calcul de Smith-Waterman).
double maximum_K (int i, int j, int n, int m, double** mat); // renvoie le max parmi la ligne gauche de mat[i][j], en fonction d'un certain coefficient de d�l�tion.
double maximum_L (int i, int j, int n, int m, double** mat); // renvoie le max parmi la colonne du haut de mat[i][j], en fonction d'un certain coefficient de d�l�tion.
void FillMatrix (int n, int m, char s1[n], char s2[m], double** mat); // remplit la matrice de similarit�.
void ShowMatrix (int n, int m,char s1[n], char s2[m], double** mat);  // affiche la matrice de similarit�.

/*****************************  Processus de tra�age et alignement  ********************************/


cell indicesMaxK (int i ,int j, int n, int m, double** mat); // r�cup�ration des indices (i-r,j) du max & sa valeur parmi la ligne � gauche de mat[i][j], en fonction d'un certain coefficient de d�l�tion.
cell indicesMaxL (int i, int j, int n, int m, double** mat); // r�cup�ration des indices (i,j-r) du max & sa valeur parmi la colonne en haut de mat[i][j], en fonction d'un certain coefficient de d�l�tion.
cell_max get_max(int n, int m,double** mat); // r�cup�ration du maximum et comptage de son nombre d'occurence. On stocke ces informations dans cell_max.
cell precell (int i, int j, int n, int m, char s1[n], char s2[m], double** mat); // calcul de la cellule d'origine/parent.
cell_liste traceback (cell c, int n, int m, char s1[n], char s2[m], double** mat); // renvoie le traceback.
void empiler(cell_liste *pile, cell p); // ajoute une cellule dans la pile qui contient le traceback.                                                                      
cell depiler (cell_liste *pile); // supprime une cellule dans la pile qui contient le traceback.
void alignement (cell_liste *chaine, int n, int m, char s1[n], char s2[m]); // affichage des alignements par rapport � notre pile.



/***********************************           MAIN             **************************************/

/************************************  Importation de fichiers  **************************************/

int main(int argc, char *argv[])
{
    if(argc<=1) // v�rifier qu'un nom de fichier a �t� donn�.
    {
        printf("Nom de fichier manquant.\n");
        return -1;  // sortir du programme le cas �ch�ant.
    }
    FILE* pf = fopen(argv[1],"rt"); 
    if(pf==NULL) // v�rification de existence du fichier.
    {
        printf("Erreur dans l'ouverture du fichier.\n");
        exit(1); // forcer la sortie.
    }

    char c;
    int l1=0,l2=0; // on compte le nombre de caract�res des cha�nes S1 et S2.
    while((c=fgetc(pf))!='\n' && c!=EOF)
        l1++;
    while((c=fgetc(pf))!='\n' && c!=EOF)
        l2++;

    rewind(pf); // Le curseur retourne au d�but du fichier.

    //on remplit nos 2 cha�nes en relisant le fichier.
    char s1[l1+1];
    int i=0;
    while((c=fgetc(pf))!='\n' && c!=EOF)
    {
        s1[i]=c; // insertion des caract�res dans la chaine s1.
        i++;
    }
    s1[l1]='\0'; // ne pas oublier le caract�re de fin de chaine.

    char s2[l2+1];
    i=0;
    while((c=fgetc(pf))!='\n' && c!=EOF)
    {
        s2[i]=c; // insertion des caract�res dans la chaine s2.
        i++;
    }
    s2[l2]='\0'; // ne pas oublier le caract�re de fin de chaine.


    // on v�rifie qu'aucune cha�ne n'est vide
    if(strcmp(s1,"")==0 || strcmp(s2,"")==0)
    {
        printf("Une chaine est vide.\n");
        exit(1); // on sort du programme.
    }

    printf("\n");
    printf("Input :");
    printf("\nSequence 1 (%d) : %s\nSequence 2 (%d) : %s\n", l1,s1, l2,s2); //on affiche nos 2 cha�nes.

    fclose(pf); //fermeture du fichier. 

/************************************ Calcul de la matrice de notation  **************************************/

    double ** mat;
    mat = allocation(l1+1,l2+1); // allocation dynamique de notre matrice de notation (l1+1 lignes et l2+1 colonnes).

    // on initialise la premi�re ligne et la premi�re colonne � 0.
    for(int i=0; i<l1+1; i++)
    {
        mat[i][0]=0.0;
    }
    for(int j=0; j<l2+1; j++)
    {
        mat[0][j]=0.0;
    }
    FillMatrix(l1+1,l2+1,s1,s2,mat); //on remplit la matrice.
    printf("\nMatrice de notation - mat[i][j] : \n\n");
    ShowMatrix(l1+1,l2+1,s1,s2,mat); //on affiche la matrice.

    // on cherche la valeur maximale de la matrice et son nombre d'occurences.
    cell_max mc = get_max(l1+1,l2+1,mat);
    printf("\nLe degre de similarite maximal est  %.1f . Nombre occurences = %d.\n",mc.val_max,mc.count);
    printf("\n");

/********************************** Calcul du traceback & affichage des alignements ******************************/
    
    /* parcours de la matrice de notation � la recherche des diff�rentes
    occurences de la valeur maximale afin d'enclencher le processus
    de traceback */

	cell cMax;
	for (int i=1; i<l1+1; i++)
		for (int j=1; j<l2+1; j++)
			if (mat[i][j] == mc.val_max) // mc.val_max valeur maximale de la matrice stock�e dans cell_max.
			{
				cMax.idx = i; cMax.idy = j; // stockage des indices dans la cell cMax.
                //printf("\nTraceback : \n");
				cell_liste chaine = traceback(cMax,l1+1,l2+1,s1,s2,mat); // pour chaque occurence du max on empile les coordon�es du traceback.
				printf("\n");
				alignement(&chaine,l1+1,l2+1,s1,s2); // on d�pile les coordonn�es du chemin pour l'affichage des alignements.
			}

    // Lib�ration de la m�moire de notre matrice (tableau � deux dimension 2).
    // Lib�ration un � un de tous les tableaux avec une boucle (l2 colonnes) puis lib�rer le tableau de pointeurs (l1 lignes).
    for (int i=0; i<l1+1; i++) {
        if(l2>0)
            free(mat[i]);
    }
    if(l1>0)
        free(mat);

    return 0;
}


/***********************************          FIN MAIN             **************************************/


/***********************************           CODE             **************************************/


/*************************************** Matrice de notation ***************************************/

/* Un tableau de double est de type double* et un tableau de tableaux de double est de type (double *)*   */

double** allocation(int l1, int l2) // allocation dynamique de notre matrice de notation. 
{
    double ** mat; // d�claration d'un tableau de dimension 2 (notre matrice).

    /* Pour allouer le tableau de l1 lignes et l2 colonnes, on commence par allouer un tableau de l1 pointeurs */
    /* on alloue ensuite chacun de ces pointeurs avec l2 double dans une boucle */

    mat = (double**) malloc ((l1+1)*sizeof(double*)); // allocation d'un tableau de l1 pointeurs. Les lignes de notre matrice � deux dimensions.
    if (mat == NULL) // v�rification de existence de la matrice.
    {
        printf("Erreur dans la matrice.\n");
        exit(1);
    }
        for (int i=0; i<l1+1; i++) {
            mat[i] = (double*) malloc ((l2+1)*sizeof(double)); // allocation d'un tableau de l2 double au bout de chaque pointeur. Pour les l1 lignes on cr�er un tableau de l2 colonnes. 
          if (mat[i] == NULL)                                               
          {
              printf("Erreur dans la matrice.\n");
              exit(1);
          }
        }
    return mat;
}

double similarite (char s1, char s2) // renvoie un coefficient de similarit� entre deux caract�res. 
{
    if (s1==s2)
        return 1.0;
    return -0.33;
}

double maxQuatrecoeff (double a, double b, double c, double d) // renvoie le coefficient maximal entre 4 valeurs (selon la m�thode de calcul de Smith-Waterman).
{
    if(a>=b && a>=c && a>=d)
    {
        return a;
    }
    if(b>=a && b>=c && b>=d)
    {
        return b;
    }
    if(c>=a && c>=b && c>=d)
    {
        return c;
    }
    return d;
}

double maximum_K (int i, int j, int n, int m, double** mat) // renvoie le max parmi la ligne gauche de mat[i][j], en fonction d'un certain coefficient de d�l�tion.
{
    double val, max = 0;
    for (int k=1; k<i; k++)
    {
        val=mat[i-k][j] - 1 - 0.33*k;
        if(max<val)
            max = val;
    }
    return max;
}

double maximum_L (int i, int j, int n, int m, double** mat) // renvoie le max parmi la colonne du haut de mat[i][j], en fonction d'un certain coefficient de d�l�tion.
{
    double val, max = 0;
    for (int l=1; l<j; l++)
    {
        val=mat[i][j-l] - 1 - 0.33*l;
        if(max<val)
            max = val;
    }
    return max;
}

void FillMatrix (int n, int m, char s1[n], char s2[m], double** mat) // remplit la matrice de similarit�.
{
    double a,b,c,d=0;
    for (int i=1; i<n+1; i++)
    {
        for (int j=1; j<m+1; j++)
        {
            a = mat[i-1][j-1] + similarite (s1[i-1],s2[j-1]); //score de l'alignement ai et bj (% � la diagonale)
            b = maximum_K(i,j,n,m,mat); //score ligne gauche
            c = maximum_L(i,j,n,m,mat); //score colonne du haut
            mat[i][j] = maxQuatrecoeff(a,b,c,d); //r�cup�ration du score maximum parmi ces 4 valeurs.
        }
    }
}

void ShowMatrix (int n, int m, char s1[n], char s2[m], double** mat) // afficher la matrice de similarit�.
{
    printf("       ");
    for (int j=1; j<m; j++) printf("%c   ",s2[j-1]); //on affiche la chaine S2
		printf("\n");

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
					if(j==0 && i>0)
						printf("%c ",s1[i-1]); // on affiche la chaine S1 � chaque d�but de ligne
					else if(j==0 && i==0)
						printf("  ");

		      printf("%.1f ", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

cell_max get_max(int n , int m, double** mat) // r�cup�ration du maximum et comptage de son nombre d'occurence. On stocke ces informations dans cell_max.
{
    cell_max mc;
    int i, j;
    mc.count= 0; 
    double max=mat[0][0];

    for(i=1; i<n; i++)
    {
        for(j=1; j<m; j++)
        {
            if(mat[i][j]>max)
            {
                max = mat[i][j];
            }
        }
    }
    mc.val_max = max;


    for(i=1; i<n; i++)
    {
        for(j=1; j<m; j++)
        {
            if(mat[i][j]==mc.val_max)
            {
                mc.count++;
            }
        }
    }

    return mc;
}

/********************************************* Processus de tra�age et alignement ******************************************/

cell indicesMaxK (int i,int j, int n,int m, double** mat) // r�cup�ration des indices (i-r,j) du max & sa valeur parmi la ligne � gauche de mat[i][j], en fonction d'un certain coefficient de d�l�tion.
{
    cell indices_maxK;
    int r=1;
    double max=0;
    double val;

    for (int k=1; k<i; k++)
    {
        val = mat[i-k][j] - 1 - 0.33*k;
        if (max<val)
        {
            max = val;
            r = k;
        }
    }

    indices_maxK.idx = i-r;
    indices_maxK.idy = j;
    indices_maxK.value = max;

    return indices_maxK;
}

cell indicesMaxL (int i,int j, int n,int m, double** mat)  // r�cup�ration des indices (i,j-r) du max & sa valeur parmi la colonne haute de mat[i][j], en fonction d'un certain coefficient de d�l�tion.
{
    cell indices_maxL;
    int r=1;
    double max=0;
    double val;
    for (int l=1; l<j; l++)
    {
        val = mat[i][j-l] - 1 - 0.33*l;
        if (max<val)
        {
            max = val;
            r = l;
        }
    }

    indices_maxL.idx = i;
    indices_maxL.idy = j-r;
    indices_maxL.value = max;

    return indices_maxL;
}

cell precell (int i,int j, int n,int m, char s1[n], char s2[m], double** mat) // calcul de la cellule d'origine/parent.
{
    cell origine;

    // recalculer les valeurs qui ont amen�es � mat[i][j].
    double a = mat[i-1][j-1] + similarite(s1[i-1], s2[j-1]); // similarit� par rapport � la diagonale.
    double d = 0;

    cell bb = indicesMaxK(i,j, n,m,mat); // r�cup�ration des indices (i-r,j) du maximum de la ligne.
    cell cc = indicesMaxL(i,j, n,m,mat); // r�cup�ration des indices (i,j-r) du maximum de la colonne.

    double prevscore = maxQuatrecoeff(a, bb.value, cc.value, d); 

    if (prevscore == a) // l'origine est � la diagonale.
    {
        origine.idx = i-1; 
        origine.idy = j-1;
        origine.value = mat[origine.idx][origine.idy];
    }
    else if (prevscore == bb.value) // origine se situe dans la ligne � gauche. 
    {
        // position � gauche avec un pas de 1.
        origine.idx = bb.idx; // p.idx = i-r;
        origine.idy = bb.idy; // p.idy = j;
        origine.value = mat[origine.idx][origine.idy];
    }
    else if (prevscore == cc.value) // origine se situe au dessus dans la colonne.
    {
        // position en haut avec un pas de 1.
        origine.idx = cc.idx;  // p.idx = i;
        origine.idy = cc.idy;  // p.idy = j-r;
        origine.value = mat[origine.idx][origine.idy];
    }
    else //cas o� mon origine vient du 0 (ce qui n'arrive jamais).
    {
        origine.idx = 0;
        origine.idy = 0;
        origine.value = 0;
    }
	
    return origine;
}

void empiler(cell_liste *pile, cell p) // ajoute une cellule dans la pile qui contient le traceback.
{
    path* c = (path*) malloc(sizeof(path)); // cr�ation d'un nouvel �l�ment de la cell_liste chain�e
    if(c==NULL) //v�rification de l'allocation en m�moire
    {
        printf("Erreur dans l'allocation de la pile !\n");
        exit(1);   
    }
    c->e = p; // ajout de l'�l�ment � empiler.
    c->next = pile->tete; // insertion en t�te de cell_liste.
    pile->tete = c; // mise � jour de la t�te de cell_liste.
}

cell_liste traceback (cell c, int n,int m, char s1[n], char s2[m], double** mat) // renvoie le traceback.
{
    cell_liste chaine;
    chaine.tete = NULL; // initalisation de la pile � vide.

    while (c.value != 0) // condition de sortie : similarit� nulle.
    {
        empiler(&chaine, c); // on ajoute cette cellule dans la pile.
        //printf("(idx , idy) = ( %d,%d ) ; Valeur =  %.1f\n", c.idx,c.idy, mat[c.idx][c.idy]);
        c = precell(c.idx,c.idy, n,m,s1,s2,mat); // on calcul le pr�decesseur de notre cellule (origine).
    }

    return chaine;
}

cell depiler (cell_liste *pile) // supprime une cellule dans la pile qui contient le traceback.
{

    path* supp_cell = pile->tete; // on r�cup�re la t�te de la cell_liste. M�morisation de l'adresse.
    if(pile == NULL || pile->tete ==NULL)
    {
        printf("Erreur au niveau de la pile.\n");
        exit(1); // forcer la sortie du programme;
    }
    cell p = supp_cell->e; // on r�cuperer les informations de la t�te de la pile.
    pile->tete = supp_cell->next; // passage au suivant avant destruction.
    free(supp_cell); // destruction de la cell m�moris�e.

    return p; // renvoie la cell en t�te de cell_liste
}

void alignement(cell_liste *chaine,int n, int m, char s1[n],char s2[m]) // affichage des alignements par rapport � notre pile.
{

    char c1[n+m]; // n+m longueur maximal possible pour un alignement dans le pire des cas.
    char c2[n+m]; //cha�nes de caract�res c1 et c2 pour nos 2 alignements.
	
    strcpy(c1 , "");  strcpy(c2 , "");  //on initialise � vide nos cha�nes. 

    int d1 = chaine->tete->e.idx; //debut de la comparaison
    int d2 = chaine->tete->e.idy; // r�cup�ration des premi�res coordon�es de la pile.

    while(chaine->tete!=NULL) // tant que la pile n'est pas vide
    {
        cell d = depiler(chaine); // on d�pile
        // d.idx et d.idy r�cup�rent les coordonn�es de la cellule qui a �t� retir�.
        while(d1<d.idx) //si il y a une d�l�tion sur s1 avant l'alignement (d.idx,d.idy)
        {
            char a[2] = {s1[d1-1]};
            strcat(c1,a);
            strcat(c2,"-");
            d1++;
        }

        while(d2<d.idy) //si il y a une d�l�tion sur s2 avant l'alignement (d.idx,d.idy)
        {
            char b[2] = {s2[d2-1]};
            strcat(c2, b);
            strcat(c1,"-");
            d2++;
        }

        if(d1==d.idx && d2==d.idy) //alignement (d.idx,d.idy)
        {
            char a[2] = {s1[d1-1]};
            char b[2] = {s2[d2-1]};
            strcat(c1,a);
            strcat(c2,b);
            d1++; d2++;
        }
    }

    printf("Output :\n");
    printf("%s\n%s\n", c1,c2); //on affiche les alignements.
}


