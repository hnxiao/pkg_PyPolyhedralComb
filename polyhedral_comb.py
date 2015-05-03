def fractional_matching_polytope(g):
# INPUT:
#     g - a graph;
# OUTPUT:
#     p - the fractional matching polytope defined on graph g.
    lp = MixedIntegerLinearProgram() # solver="PPL"
    x = lp.new_variable(real=True, nonnegative=True)
    for v in g:
        edges = g.edges_incident(v, labels = False)
        lp.add_constraint( sum([x[Set(e)] for e in edges]) <= 1 ) # matching constraint
    ps = list(Subsets(Set(range(len(g.vertices()))))) # powerset of vertex set V
    for U in ps:
        if U.cardinality()%2 ==1: # odd subset of vertex set V
            edges = g.subgraph(U).edges(labels = False)
            if edges !=[]:
                lp.add_constraint( sum([x[Set(e)] for e in edges]) <= (U.cardinality()-1)/2 ) # Edmonds odd set constraint
    p=lp.polyhedron()
    return p


def preference_list(g):
# INPUT:
#     g - a graph;
# OUTPUT:
#     pl - a random preference list on vertices of g.
    pl=[];
    for v in g:
        pl.append(Permutations(g.neighbors(v)).random_element())
    return pl

def fractional_stable_matching_polytope(g, pl):
# Fractional stable matching polytope defined by Rothblum system
# INPUT:
#     g - a graph;
#     pl - a preference list of graph g.
# OUTPUT:
#     P - the fractional stable matching polytope defined on g with preference list pl.
    lp = MixedIntegerLinearProgram() # solver="PPL"
    x = lp.new_variable(real=True, nonnegative=True)
    for v in g: # Incidence Matrix
        edges = g.edges_incident(v, labels = False)
        lp.add_constraint( sum([x[Set(e)] for e in edges]) <= 1 ) # matching constraint
    for (u,v) in g.edges(labels=False):
        lp.add_constraint(
            sum(x[Set((u,i))] for i in g.neighbors(u) if pl[u].index(i) > pl[u].index(v))
            + sum(x[Set((v,j))] for j in g.neighbors(v) if pl[v].index(j) > pl[v].index(u))
            + x[Set((u,v))] >= 1) # Rothblum stability constraint
    p=lp.polyhedron()
    return p


def fractional_stable_matching_polytope_RE(g, pl):
# Fractional stable matching polytope defined by Rothblum-Edmonds system
# INPUT:
#     g - a graph;
#     pl - a preference list of graph g.
# OUTPUT:
#     P - the fractional stable matching polytope defined on g with preference list pl.
    lp = MixedIntegerLinearProgram() # solver="PPL"
    x = lp.new_variable(real=True, nonnegative=True)
    for v in g: # Incidence Matrix
        edges = g.edges_incident(v, labels = False)
        lp.add_constraint( sum([x[Set(e)] for e in edges]) <= 1 ) # matching constraint
    ps = list(Subsets(Set(range(len(g.vertices()))))) # powerset of vertex set V
    for U in ps:
        if U.cardinality()%2 ==1: # odd subset of vertex set V
            edges = g.subgraph(U).edges(labels = False)
            if edges !=[]:
                lp.add_constraint( sum([x[Set(e)] for e in edges]) <= (U.cardinality()-1)/2 ) # Edmonds odd set constraint
    for (u,v) in g.edges(labels=False):
        lp.add_constraint(
            sum(x[Set((u,i))] for i in g.neighbors(u) if pl[u].index(i) > pl[u].index(v))
            + sum(x[Set((v,j))] for j in g.neighbors(v) if pl[v].index(j) > pl[v].index(u))
            + x[Set((u,v))] >= 1) # Rothblum stability constraint
    p=lp.polyhedron()
    return p


def fractional_feedback_vertex_set_polyhedron(d):
# INPUT:
#    d - a digraph;
# OUTPUT:
#    p - the fractional feedback vertex set polytope defined on digraph d
    lp = MixedIntegerLinearProgram() # solver="PPL"
    x = lp.new_variable(real=True, nonnegative=True)
    cs = d.all_simple_cycles() # set of all simple cycles
    for i in range(len(cs)):
        lp.add_constraint( sum([x[v] for v in cs[i]]) >= 1 )   # cycle-vertex matrix
    p = lp.polyhedron()
    return p


def fractional_feedback_arc_set_polyhedron(d):
# Another implementation by testing each subset S of arc set A and check whethr G-S is acyclic?
# INPUT:
#    d - a directed graph;
# OUTPUT:
#    p - the fractional feedback arc set polytope defined on digraph d
    lp = MixedIntegerLinearProgram() # solver="PPL"
    x = lp.new_variable(real=True, nonnegative=True)
    cs = d.all_simple_cycles() # set of all simple cycles
    for i in range(len(cs)):
        lp.add_constraint( sum([x[a] for a in d.subgraph(cs[i]).hamiltonian_cycle().edges(labels=False)]) >= 1 )   # cycle-arc matrix
    p = lp.polyhedron()
    return p


def feedback_vertex_sets_minimal(d):
# This function depends on functin feedback_vertex_set_polyhedron and returns all the minimal feedback vertex sets
# INPUT:
#    d - a digraph;
# OUTPUT:
#    s - set of all minimal feedback vertex sets of digraph d
    s = []
    vl = feedback_vertex_set_polyhedron(d).vertices_list()
    for i in range(len(vl)):
        if all(e in ZZ for e in vl[i]): # all integral vertices of feedback vertex set polyhedron
            s.append(vl[i])
    return s


def fractional_independent_set_polyhedron(g):
# INPUT:
#    g - a graph;
# OUTPUT:
#    p - the fractional independent set polyhedron
    lp = MixedIntegerLinearProgram()
    x = lp.new_variable(real=True, nonnegative=True)
    for (u,v) in g.edges(labels = False):
        lp.add_constraint( x[u] + x[v] <= 1 ) # independence constraint
    p=lp.polyhedron()
    return p

def fractional_dominating_set_polyhedron(g):
# INPUT:
#    g - a graph;
# OUTPUT:
#    p - the fractional dominating set polyhedron
    lp = MixedIntegerLinearProgram()
    x = lp.new_variable(real=True, nonnegative=True)
    for u in g:
        p.add_constraint( x[u] + sum([x[v] for v in g.neighbors(u)]) >= 1 ) # domination constraint
    p = lp.polyhedron()
    return p


def fractional_kernel_polytope(d):
# A generic function that generates fractional kernel polytopes
# INPUT:
#    d - an acyclic digraph;
# OUTPUT:
#    p - the fractional kernel polytope defined on digraph d
    lp = MixedIntegerLinearProgram()
    x = lp.new_variable(real=True, nonnegative=True)
    for u in d:
        lp.add_constraint( x[u] + sum([x[v] for v in d.neighbors_out(u)]) >= 1 ) # domination constraint
    if d.is_directed_acyclic():
        for (u,v) in d.edges(labels = False):
            lp.add_constraint( x[u] + x[v] <= 1 ) # acyclic independence constraint
    else:
        cs = d.to_undirected().cliques_maximal() # set of all maximal cliques
        for clique in cs:
            lp.add_constraint( sum(x[v] for v in clique) <= 1 ) # general independence constraint
    p = lp.polyhedron()
    return p

"""
def fractional_kernel_polytope_acyclic(d): # ayclic definition is from Egres Open
# INPUT:
#    d - an acyclic digraph;
# OUTPUT:
#    p - the fractional kernel polytope defined on acyclic digraph d
    lp = MixedIntegerLinearProgram()
    x = lp.new_variable(real=True, nonnegative=True)
    for (u,v) in d.edges(labels = False):
        lp.add_constraint( x[u] + x[v] <= 1 )
    for u in d:
        lp.add_constraint( x[u] + sum([x[v] for v in d.neighbors_out(u)]) >= 1 )
    p = lp.polyhedron()
    return p

def fractional_kernel_polytope(d):
# INPUT:
#    d - an acyclic digraph;
# OUTPUT:
#    p - the fractional kernel polytope defined on digraph d
    lp = MixedIntegerLinearProgram()
    x = lp.new_variable(real=True, nonnegative=True)
    cs = d.to_undirected().cliques_maximal() # the set of all maximal cliques
    for clique in cs:
        lp.add_constraint( sum(x[v] for v in clique) <= 1 )
    for u in d:
        lp.add_constraint( x[u] + sum([x[v] for v in d.neighbors_out(u)]) >= 1 )
    p = lp.polyhedron()
    return p
"""