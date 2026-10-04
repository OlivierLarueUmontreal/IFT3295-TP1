"""IFT3295 - TP1 - Assemblage de fragments (Ex. 2.2 et 2.3).

Ce fichier est A COMPLETER. Regles:
* Completez uniquement le corps des fonctions ci-dessous.
* Ne modifiez ni les noms, ni les signatures, ni les types.
* Vous pouvez utiliser NetworkX pour representer et manipuler le graphe.
* Typez tous vos parametres et retours (types simples: int, str, list, ...).
* Le programme principal est ``main.py`` (fourni): il appelle les fonctions
  de ce module, il ne faut donc rien executer ici directement.
"""

import networkx as nx

SEUIL_DEFAUT: int = 80

# Les noeuds du graphe sont les entiers 0..n-1 (indices des reads); chaque
# arete porte l'attribut entier ``poids`` (= score du chevauchement).


def graphe_chevauchements(
    scores: list[list[int]], seuil: int = SEUIL_DEFAUT
) -> nx.DiGraph:
    """Construit le graphe oriente des chevauchements pertinents.

    Args:
        scores (list[list[int]]): Matrice carree des scores de chevauchement.
        seuil (int): Score minimal requis pour conserver une arete.

    Returns:
        nx.DiGraph: Graphe dont les noeuds sont les indices ``0..n-1``.
        Une arete ``(i, j)`` existe si ``scores[i][j] >= seuil``; son attribut
        ``poids`` vaut le score correspondant.
    """
    raise NotImplementedError  # TODO


def reduction_transitive(graphe: nx.DiGraph) -> nx.DiGraph:
    """Calcule la reduction transitive du graphe de chevauchement.

    Supprime une arete ``(u, v)`` lorsqu'un autre chemin simple de ``u`` a
    ``v`` subsiste. Pour un cycle de deux noeuds, conserve uniquement l'arete
    de plus grand poids.

    Args:
        graphe (nx.DiGraph): Graphe oriente dont les aretes portent un attribut
            entier ``poids``.

    Returns:
        nx.DiGraph: Copie reduite du graphe d'entree. Le graphe d'entree n'est
        pas modifie.
    """
    raise NotImplementedError  # TODO
def ordre_assemblage(graphe: nx.DiGraph) -> list[int]:
    """Trouve l'ordre d'assemblage des reads dans le graphe reduit.

    Args:
        graphe (nx.DiGraph): Graphe reduit produit par
            ``reduction_transitive``.

    Returns:
        list[int]: Indices des reads dans un chemin du graphe qui visite
        chaque noeud exactement une fois. L'enonce garantit l'existence de
        ce chemin apres reduction.
    """
    raise NotImplementedError  # TODO


def sequence_finale(reads: list[str], ordre: list[int]) -> tuple[str, list[int]]:
    """Assemble les reads en une sequence de fragment genomique.

    Args:
        reads (list[str]): Sequences des reads a assembler.
        ordre (list[int]): Indices des reads dans l'ordre d'assemblage.

    Returns:
        tuple[str, list[int]]: Tuple contenant la sequence assemblee et les
        longueurs des chevauchements entre les paires de reads consecutifs,
        dans l'ordre.
    """
    raise NotImplementedError  # TODO
