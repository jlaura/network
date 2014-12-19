import numpy as np

def networkg(ntw, pointpattern, nsteps=10, permutations=99, threshold=5, envelope=False, poisson=False)
    if not pointpattern in ntw.pointpatterns.keys():
        print "Key Error: Available point patterns are {}".format(ntw.pointpatterns.keys())
        return


    if poisson:
        distribution = 'poisson'
    else:
        distribution = 'uniform'
    npts = len(nts.pointpatterns[pointpattern])

    results = np.empty((npts, permutations + 1))

    for p in xrange(permutations):
        rpts = ntw.simulate_observations(ntw,npts, distribution=distribution)

        dist_to_node = {}
        c = 0
        for edge, coords in rpts.iteritems():
            for p in coords:
                #Compute the distance from the new point to the nodes
                d1, d2 = compute_distance_to_nodes(ntw, p[0], p[1], edge)
                dist_to_node[c] = {edge[0]:d1, edge[1]:d2}
                c += 1

        #2 column matrix, c0 is the nn id and c1 is the distance

        #Move to util.py
        nearest = np.zeros((len(pt_indices), 2), dtype=np.float32)

        for i, p1 in enumerate(pt_indices):
            dist1, dist2 = dist_to_node[p1].values()
            endnode1, endnode2 = dist_to_node[p1].keys()

            snapped_coords = pts.snapped_coordinates[p1]
            nearest_obs1, nearest_node1, nearest_node_distance1 = nearestneighborsearch(obs_to_node, alldistances, endnode1, dist1)
            nearest_obs2, nearest_node2, nearest_node_distance2 = nearestneighborsearch(obs_to_node, alldistances, endnode2, dist2)

            if nearest_node_distance2 <= nearest_node_distance1:
                nearest[i,0] = nearest_obs2
                nearest[i,1] = nearest_node_distance2
            else:
                nearest[i,0] = nearest_obs1
                nearest[i,1] = nearest_node_distance1

        xsim, ysim = gstat(nearest, nsteps=100)



def gfunction(nearest):
    """
    Compute a G-Function

    Parameters
    ----------
    nsteps          int The number of distance bands
    permutations    int The number of permutations to perform
    threshold       int Upper and lower significance threshold
    envelope        bool Return results of all permutations
    poisson         bool Use a poisson distribution to
                         determine the number of points
    """
    lowerbound = np.min(nearest[:,1])
    upperbound = np.max(nearest[:,1])
    nobs = len(nearest[:,1])
    x = np.linspace(lowerbound, upperbound, nsteps)
    nearest = np.sort(nearest[:,1])

    y = np.empty(len(x))
    for i,r in enumerate(x):
        g = len(np.where(nearest <= r)[0]) / float(nobs)
        y[i].append(g)

    return x, y

    npts = len(nearest[:,0])
    for p in permutations:
        npts = len(

