from itertools import chain


def possible_ltracks(lcts):
    MEp = ([a, b, c, d] for a in lcts[(1, 1)] for b in lcts[(1, 2)]
           for c in lcts[(1, 3)] for d in lcts[(1, 4)])
    MEm = ([a, b, c, d] for a in lcts[(-1, 1)] for b in lcts[(-1, 2)]
           for c in lcts[(-1, 3)] for d in lcts[(-1, 4)])
    return chain(MEp, MEm)


def ltrack_orientation(event, track):
    pass


STATION_Z = {
    (-1, 4): -10.5,
    (-1, 3): -9.5,
    (-1, 2): -8,
    (-1, 1): -6,
    (1, 1): 6,
    (1, 2): 8,
    (1, 3): 9.5,
    (1, 4): 10.5
}


def lct_xyz(event, lct):
    theta = event.hit_theta[lct]
    phi = event.hit_phi[lct]
    z = STATION_X[(event.hit_endcap[lct], event.hit_station[lct])]
    y = tan(theta) * cos(phi) * z
    x = tan(theta) * cos(phi) * z
