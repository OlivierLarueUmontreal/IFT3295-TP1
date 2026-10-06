"""IFT3295 - TP1 - Algorithme de chevauchement de sequences (Q1.5, Ex. 2.1).

Ce fichier est A COMPLETER. Regles:
* Completez uniquement le corps des fonctions ci-dessous.
* Ne modifiez ni les noms, ni les signatures, ni les types, ni les constantes.
* Typez tous vos parametres et retours (types simples: int, str, list, ...).
* Le programme principal est ``main.py`` (fourni): il appelle les fonctions
  de ce module, il ne faut donc rien executer ici directement.
"""

MATCH: int = 4
MISMATCH: int = -4
INDEL: int = -8


def chevauchement_maximal(x: str, y: str) -> tuple[int, str, str, int]:
    """Calcule le meilleur chevauchement ordonne entre deux sequences.

    Args:
        x (str): Premiere sequence a aligner.
        y (str): Deuxieme sequence a aligner.

    Returns:
        tuple[int, str, str, int]: Score maximal, deux lignes alignees et
        longueur du chevauchement.

    Examples:
        >>> chevauchement_maximal("ACCA", "CACGC")
        (8, 'CA', 'CA', 2)
        >>> chevauchement_maximal("CACGC", "ACCA")
        (4, 'ACGC', 'AC-C', 4)
    """
    n: int = len(x)
    m: int = len(y)

    # + 1 pour les case 0 de la matrice
    V: list[list[int]] = [[0] * (m + 1) for _ in range(n + 1)]

    # Puisque l'on a toujours un prefix de y
    for j in range(1, m + 1):
        V[0][j] = j * INDEL

    # V[i][0] = 0 pour tout i de 1 a n + 1 puisque qu'on prend suffix de x, donc i premiers charactere ne sont jamais pris

    # Meme recurrence que pour l'alignement global, O(nm)
    # basé sur: https://fr.wikipedia.org/wiki/Algorithme_de_Needleman-Wunsch
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            delta: int = MATCH if x[i - 1] == y[j - 1] else MISMATCH
            diagonal: int = V[i - 1][j - 1] + delta
            left: int = V[i - 1][j] + INDEL
            up: int = V[i][j - 1] + INDEL
            V[i][j] = max(diagonal, up, left)

    # Cherche le score max pour backtrack, V[n][j] pour s'assurer d'etre a la fin de x
    best_j: int = 0
    for j in range(1, m + 1):
        if V[n][j] > V[n][best_j]:
            best_j = j
    score: int = V[n][best_j]

    # backtracking
    alignement_x: list[str] = []
    alignement_y: list[str] = []
    i: int = n
    j: int = best_j

    while i > 0 and j > 0:
        current = V[i][j];
        diagonal = V[i - 1][j - 1]
        left = V[i-1][j]
        up = V[i][j-1]
        delta = MATCH if x[i - 1] == y[j - 1] else MISMATCH

        if current == diagonal + delta:
            alignement_x.append(x[i-1])
            alignement_y.append(y[j-1])
            i -= 1
            j -= 1
        elif current == left + INDEL:
            alignement_x.append(x[i-1])
            alignement_y.append("-")
            i -= 1
        elif current == up + INDEL:
            alignement_x.append("-")
            alignement_y.append(y[j-1])
            j -= 1

    # Comme i atteint 0 avant j (puisque suffix de x et prefix de y) on ajoute des - pour la fin des j (x serait des - aligné avec un char de y)
    while j > 0:
        alignement_x.append("-")
        alignement_y.append(y[j-1])
        j-=1

    # Resverse l'alignement et retourne en string
    ligne_x: str = "".join(reversed(alignement_x))
    ligne_y: str = "".join(reversed(alignement_y))
    return score, ligne_x, ligne_y, len(ligne_x)


def matrice_chevauchements(reads: list[str]) -> list[list[int]]:
    """Construit la matrice des scores de chevauchement de toutes les paires.

    Args:
        reads (list[str]): Sequences a comparer, par exemple les 20 reads de
            ``reads.fq``.

    Returns:
        list[list[int]]: Matrice carree dont l'element ``[i][j]`` est le score
        maximal de la paire ordonnee ``(reads[i], reads[j])``. La diagonale
        contient des zeros.
    """
    n: int = len(reads)
    scores: list[list[int]] = [[0] * n for _ in range(n)]

    # calcul le score max pour chaque pair i,j
    for i in range(n):
        for j in range(n):
            if i != j: # Pour éviter la diagonale nulle
                scores[i][j] = chevauchement_maximal(reads[i], reads[j])[0]

    return scores
