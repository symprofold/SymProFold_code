

def permutations(models0, models1):
    """ Generate all possible permutations of models0 and models1. """

    permutations = []

    for i in range(0, max([max(models0), max(models1)])+1):

        for j,m0 in enumerate(models0):
            if m0 > i:
                continue
            
            for k,m1 in enumerate(models1):
                if m1 > i:
                    continue

                permutation = (m0, m1)

                if permutation not in permutations:
                    permutations.append(permutation)

    return permutations
