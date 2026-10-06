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

from overlap import chevauchement_maximal

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
    # Basé sur la doc NetworkX: https://networkx.org/documentation/stable/reference/classes/digraph.html
    graphe = nx.DiGraph()
    n = len(scores)

    graphe.add_nodes_from(range(n))

    for i in range(n):
        for j in range(n):
            if i != j and scores[i][j] >= seuil:
                graphe.add_edge(i, j, poids=scores[i][j])

    return graphe

# Fonction utilitaire me permettant de voir si un chemin existe du noeud départ a arriver existe
def chemin_existe(graphe: nx.DiGraph, depart: int, arriver: int) -> bool:
    a_visiter: list[int] = [depart]
    visiter: set[int] = {depart}

    while a_visiter:
        noeud: int = a_visiter.pop()
        if noeud == arriver:
            return True
        
        # ajoute tout les voisins du noeuds courant dans la liste a visiter
        for voisin in graphe.successors(noeud):
            if voisin not in visiter:
                visiter.add(voisin)
                a_visiter.append(voisin)

    return False

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

    graphe_reduit = graphe.copy()

    # Pour cycle de deux noeuds. Si u -> v == v -> u, je garde u -> v
    for u, v in list(graphe_reduit.edges):
        if u < v and graphe_reduit.has_edge(v, u):
            if graphe_reduit[u][v]["poids"] >= graphe_reduit[v][u]["poids"]:
                graphe_reduit.remove_edge(v, u)
            else:
                graphe_reduit.remove_edge(u, v)


    for u, v in list(graphe_reduit.edges):
        poids: int = graphe_reduit[u][v]["poids"]
        graphe_reduit.remove_edge(u, v)

        # On retire le edge initialement, s'il n'y a pas d'autre chemin possible u -> v,
        #  alors on remet l'arrete, sinon on la laisse retirer, elle est redondante
        if not chemin_existe(graphe_reduit, u, v):
            graphe_reduit.add_edge(u, v, poids)


    return graphe_reduit


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

    # Verifie que le graphe est linéaire: 1 arret sortant de chaque node
    if any(graphe.out_degree(node) > 1 for node in graphe.nodes):
        raise ValueError("Graphe reduit non lineaire")

    # Trouver le premier node: https://stackoverflow.com/questions/61845326/with-networks-how-to-find-first-nodes-in-a-digraph
    root = [node for node, in_degree in graphe.in_degree if in_degree == 0]
    if len(root) > 1 :
        raise ValueError(f"Erreur, il n'y a pas qu'un seul node de depart, mais bien: {len(depart)}")

    current = root[0]
    result = [current]
    # set pour verifier si un node a deja ete vu lors de  l'assemblage: O(1)
    # Au lieu de O(n) si je regarde dans result
    seen: set[int] = {current} 

    while True:
        successeurs = list(graphe.sucessors(current))

        if not successeurs: break # s'il n'y a plus de successeur: terminer

        current = successeurs[0]
        if current in seen: # Si on a deja vu le node, on a un cycle
            raise ValueError(f"Cycle detecter a {current}")

        seen.add(current)
        result.append(current)
        
    return result


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