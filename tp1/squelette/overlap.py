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
    raise NotImplementedError  # TODO


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
    raise NotImplementedError  # TODO
