from sage.all import *

import argparse
import itertools as it
import os
import pprint

from scipy.stats import poisson


# Maximum number of tries to use when generating matrix A.
# Each try will fail if A is not full rank, therefore we may need more than one try.
MAX_TRIES_FOR_A_GEN = 10

# Maximum number of tries to use when generating matrix V.
# Each try will fail if V does not have different rows.
MAX_TRIES_FOR_V_GEN = 10


def gen_binary_vector_of_fixed_weight(n, w):
    support = sample(range(n), w)
    vector = zero_vector(GF(2), n)
    for i in support:
        vector[i] = 1
    return vector


class BinaryPKP:

    def __init__(self, m, n, l):
        """
        Args:
            m (int): The number of rows (rank) of public matrix A
            n (int): The number of columns of matrix A
            l (int): The number of columns of Vpub and Vsec
        """
        assert(m < n)

        self.m = m
        self.n = n
        self.l = l

    def gen_Vsec_from_A(self, A, max_tries_for_V_gen=MAX_TRIES_FOR_V_GEN):
        Aker = A.right_kernel()

        for _try_V in range(max_tries_for_V_gen):
            Vsec = matrix([Aker.random_element() for _ in range(self.l)]).transpose()

            listV = list(Vsec)

            # Ensures Vpub has no two equal rows
            if all([(listV.count(row) == 1) for row in listV]):
                return Vsec

        raise RuntimeError('Could not Vsec without repeating rows.')


    def keygen_weak(self,
                    attack_param_w,
                    attack_param_la,
                    max_tries_for_A_gen=MAX_TRIES_FOR_A_GEN,
                    max_tries_for_V_gen=MAX_TRIES_FOR_V_GEN):
        """
        This function generates a random key pair for binary PKP that can be attacked with parameters
        (w, la) = (attack_param_w, attack_param_la). That is, it returns a triple (A, Vsec, Vpub) such
        that
            A is an m by n binary matrix
            V and Vsec are two n by l binary matrices
            A*Vsec = 0
            Vsec has no two equal rows (security requirement for PKP)
            Vpub is a permutation of Vsec
            The rowspace of A contains at least attack_param_la vectors of weight w

        Notice that we require that attack_param_la is limited <= m which is enough for the attack to
        work and it is much easier to generate.

        Args:
            attack_param_w (int): Attack parameter w
            attack_param_la (int): Attack parameter la
        """

        assert(attack_param_la <= self.m)

        for _try_A in range(max_tries_for_A_gen):

            A_low_weight_rows = [gen_binary_vector_of_fixed_weight(self.n, attack_param_w)
                                 for _ in range(attack_param_la)]
            A_random_rows = [random_vector(GF(2), self.n) for _ in range(attack_param_la, self.m)]


            A = matrix(GF(2), A_low_weight_rows + A_random_rows)

            if A.rank() < self.m:
                continue

            try:
                Vsec = self.gen_Vsec_from_A(A, max_tries_for_V_gen=max_tries_for_V_gen)
                Vpub = matrix(sorted(Vsec))

                return A, Vsec, Vpub

            except RuntimeError:
                pass

        raise ValueError('Could not generate a weak key with the selected parameters.')


    def keygen(self,
               max_tries_for_A_gen=MAX_TRIES_FOR_A_GEN,
               max_tries_for_V_gen=MAX_TRIES_FOR_V_GEN):
        """
        This function generates a random key pair for binary PKP.
        That is, it returns a triple (A, Vsec, Vpub) such that
            A is an m by n binary matrix
            V and Vsec are two n by l binary matrices
            A*Vsec = 0
            Vsec has no two equal rows (security requirement for PKP)
            Vpub is a permutation of Vsec

        """

        for _try_A in range(max_tries_for_A_gen):

            A = random_matrix(GF(2), self.m, self.n)

            if A.rank() < self.m:
                continue

            try:
                Vsec = self.gen_Vsec_from_A(A, max_tries_for_V_gen=max_tries_for_V_gen)
                Vpub = matrix(sorted(Vsec))

                return A, Vsec, Vpub

            except RuntimeError:
                pass

        raise ValueError('Could not generate a key with the selected parameters.')


def get_all_vectors_of_given_weight(n, w):
    vec = [0]*(n - w) + [1]*w

    while True:
        for k in range(n - 2, -1, -1):
            if vec[k] < vec[k + 1]:
                break
        else:
            # vec is the last permutation
            return

        for l in range(n - 1, -1, -1):
            if vec[k] < vec[l]:
                break

        vec[k], vec[l] = vec[l], vec[k]
        vec[k + 1:] = reversed(vec[k + 1:])

        yield(vector(GF(2), vec))


def reversed_column_permutation(x, columns):
    v = [None]*len(x)
    for i in range(len(x)):
        v[columns[i]] = x[i]

    return v


def get_G_permutation_in_systematic_form(G, max_tries=100):
    k, n = G.dimensions()

    for i in range(max_tries):

        columns = list(range(n))
        shuffle(columns)

        Gprime = G.matrix_from_columns(columns)
        Gprime.echelonize()

        if Gprime[:,:k] == identity_matrix(k):
            break

    else:
        raise RuntimeError('Could not find a permutation of G that can be written in systematic \
            form')

    return Gprime, columns


def stern_iter(G, w, p, l):
    k, n = G.dimensions()

    Gprime, cols_permutation = get_G_permutation_in_systematic_form(G)

    set_of_low_weight_vectors = set()

    H0 = {}
    for u in get_all_vectors_of_given_weight(floor(k/2), p):
        x = u*Gprime[:floor(k/2)]
        x_l = tuple(x[i] for i in range(k, k + l))

        if x_l not in H0:
            H0[x_l] = []
        H0[x_l].append(x)

    H1 = {}
    for u in get_all_vectors_of_given_weight(ceil(k/2), p):
        x = u*Gprime[-ceil(k/2):]
        x_l = tuple(x[i] for i in range(k, k + l))

        if x_l not in H1:
            H1[x_l] = []
        H1[x_l].append(x)

    for x_l in H0:
        try:
            for x0, x1 in it.product(H0[x_l], H1[x_l]):
                x = x0 + x1
                # print(x)
                if x.hamming_weight() == w:
                    v = reversed_column_permutation(x, cols_permutation)
                    set_of_low_weight_vectors.add(tuple(v))
        except KeyError:
            pass

    return set_of_low_weight_vectors


def stern(G, w, p, l, n_iters, progress=True, max_errors_systematic_form=1000, verbose=False):
    '''
    This function implements Stern's algorithm for finding low weight codewords.
    The reader will notice that it follows closely the excellent description of the algorithm given
    in Johansson and Löndahl's paper (https://lup.lub.lu.se/search/ws/files/6012849/2204899.pdf).
    '''
    import sys

    set_of_low_weight_vectors = set()

    i_iter = 0
    bulk_iter = 0
    n_errors = 0
    while(i_iter < n_iters):
        bulk_iter += 1

        try:
            found = stern_iter(G, w, p, l)
        except RuntimeError:
            n_errors += 1
            if verbose:
                print(f'Stern: error #{n_errors} while trying to systematize matrix. Continuing...')

            if n_errors >= max_errors_systematic_form:
                raise RuntimeError('Stern: reached max number of errors when trying to syst. the matrix')
            found = set()

        set_of_low_weight_vectors |= found

        i_iter += len(found)

        n_unique_codewords = len(set_of_low_weight_vectors)
        n_total_codewords = i_iter

        s = (f'[%2d%% complete] # Processed codewords: ({n_unique_codewords} unique, {n_total_codewords} total)'
                % min(100, (i_iter*100/n_iters)))

        print('\r' + s, end='')

    print()



    return [vector(v) for v in set_of_low_weight_vectors]


class Attack:

    def __init__(self, A, Vpub, attack_param_w, attack_param_la):
        self.A = A
        self.Vpub = Vpub
        self.w = attack_param_w
        self.la = attack_param_la

    def get_parameters_for_stern_algorithm(self, target_matrix):
        k, n = target_matrix.dimensions()

        def get_l(p):
            return floor(log(binomial(k//2, p)))

        def target(p):
            l = get_l(p)

            den = (binomial(n - k - l, self.w - 2*p) * binomial(k//2, p))

            if den == 0:
                return infinity

            t = 2*l*binomial(n, self.w) / den
            return numerical_approx(t)

        values = {p: target(p) for p in range(1, self.w//2 + 1)}

        p = min(values, key=lambda p: values[p])

        return (p, get_l(p))

    def get_expected_number_of_low_weight_codewords(self, target_matrix):
        prob_of_being_in_rowspace = pow(2, target_matrix.nrows() - target_matrix.ncols())
        n_low_weight_vectors = binomial(target_matrix.ncols(), self.w)

        return numerical_approx(self.la + n_low_weight_vectors * prob_of_being_in_rowspace)

    def find_low_weight_codewords_in_rowspace(self, target_matrix, error=1e-6):
        p, l = self.get_parameters_for_stern_algorithm(target_matrix)
        n_codewords = self.get_expected_number_of_low_weight_codewords(target_matrix)

        gamma = -log(error, n_codewords) + 1
        n_iters = ceil(gamma * n_codewords * log(n_codewords, 2))

        # print(f'Running Stern`s algorithm until {n_iters} (non-unique) codewords low weight codewords are processed')
        print(f'Running Stern`s algorithm with parameters p={p}, l={l}, error={error}')
        return stern(target_matrix, self.w, p, l, n_iters)

    def phase1_find_low_weight_keywords(self):
        '''
        Returns the low weight sets of vectors LWSA and LWSK, corresponding to vectors in the
        rowspace of A and K, respectively, where K is the left kernel matrix of Vpub.
        '''

        LWSA = self.find_low_weight_codewords_in_rowspace(self.A)
        K = self.Vpub.left_kernel().basis_matrix()
        LWSK = self.find_low_weight_codewords_in_rowspace(K)

        return LWSA, LWSK


def factorial_ext(x):
    if x >= 1:
        return gamma(x + 1)
    return 1

def approx_binomial(x, y):
    return factorial_ext(x)/(factorial_ext(x - y)*factorial_ext(y))

def get_column_permutation_from_B_to_A(A, B):
    '''
    Given two matrices A and B that are equivalent up to a permutation of columns, this function
    computes one permutation perm such that perm(B) = A.
    '''
    perm = []
    lA = list(A)
    for i, _ in enumerate(B):
        perm.append(lA.index(B[i]))
    return perm


class Estimates:
    '''
    This class contains the implementation of the estimates on the attack performance.
    In general, each estimate_* method corresponds to an estimation, analytical or based on
    simulations, that is described in the paper.
    '''

    def __init__(self, m, n, l):
        self.m = m
        self.n = n
        self.l = l

    def estimate_fraction_of_weak_keys(self, attack_param_w, attack_param_la):
        N = binomial(self.n, attack_param_w)
        probability_of_being_in_rowspace = 2**self.m / 2**self.n

        pois = poisson(N*probability_of_being_in_rowspace)
        return 1 - pois.cdf(attack_param_la - 1)

    def estimate_size_of_LWSK(self, attack_param_w, attack_param_la):
        av_size_of_LWSK = (
            attack_param_la + 2**(self.n - self.l) / 2**self.n * approx_binomial(self.n,
                                                                                 attack_param_w)
        )
        return numerical_approx(av_size_of_LWSK)

    def estimate_median_q_alpha_with_simulations(self, attack_param_w, alpha, nsamples=500):

        q_alphas = []
        for _ in range(nsamples):

            count = 0
            L = matrix([gen_binary_vector_of_fixed_weight(self.n, attack_param_w)
                        for _ in range(alpha)])

            Lcols_tuples = [tuple(x) for x in (L.T)]
            Lcols_w_sample = sample(Lcols_tuples, attack_param_w)

            classes = set(tuple(t) for t in Lcols_w_sample)
            counts = prod([binomial(Lcols_tuples.count(c), Lcols_w_sample.count(c))
                           for c in classes])
            prob = counts / binomial(self.n, attack_param_w)

            q_alphas.append(prob)

        return numerical_approx(median(q_alphas))

    def estimate_median_q_alpha_analitically(self, attack_param_w, alpha, nsamples=500):
        prod_k = 1
        for k in range(0, alpha + 1):
            pk = self.get_probability_p_k_alpha(attack_param_w, k, alpha)
            value = max(1,
                        pow(approx_binomial(self.n*pk, attack_param_w*pk),
                            approx_binomial(alpha, k)))
            prod_k *= numerical_approx(value)

        q_alpha = prod_k/approx_binomial(self.n, attack_param_w)
        return q_alpha

    def get_probability_p_k_alpha(self, attack_param_w, k, alpha):
        return (attack_param_w/self.n)**k * (1 - attack_param_w/self.n)**(alpha - k)

    def estimate_nchildren_per_level(self,
                                     attack_param_w,
                                     attack_param_la,
                                     with_sim=False,
                                     sim_sample=500):

        av_size_of_LWSK = self.estimate_size_of_LWSK(attack_param_w, attack_param_la)
        nchildren_per_level = [av_size_of_LWSK]

        for alpha in range(1, attack_param_la):
            if with_sim:
                q_alpha = self.estimate_median_q_alpha_with_simulations(attack_param_w,
                                                                        alpha,
                                                                        nsamples=sim_sample)
            else:
                q_alpha = self.estimate_median_q_alpha_analitically(attack_param_w, alpha)

            nchildren_per_level.append(av_size_of_LWSK*q_alpha)

        return nchildren_per_level

    def estimate_work_factor_permutations(self, attack_param_w, attack_param_la):
        perms = 1
        for k in range(0, attack_param_la + 1):
            pk = self.get_probability_p_k_alpha(attack_param_w, k, attack_param_la)
            perms *= numerical_approx(
                max(1, pow(factorial_ext(self.n*pk)/factorial_ext(self.m*pk),
                           float(binomial(attack_param_la, k)))))

        return perms

    def get_work_factor_search_given_nchildren_per_level(self, nchildren_per_level):
        return prod([x for x in nchildren_per_level if x > 1])

    def get_all_estimates_for_attack_parameters(self, attack_param_w, attack_param_la):
        data = {}
        av_size_of_LWSK = self.estimate_size_of_LWSK(attack_param_w, attack_param_la)

        nchildren_per_level_simulation = self.estimate_nchildren_per_level(attack_param_w,
                                                                           attack_param_la,
                                                                           with_sim=True,
                                                                           sim_sample=1000)

        nchildren_per_level_analytic = self.estimate_nchildren_per_level(attack_param_w,
                                                                         attack_param_la,
                                                                         with_sim=False)

        fraction_of_keys = self.estimate_fraction_of_weak_keys(attack_param_w,
                                                               attack_param_la)

        work_factor_of_search_simulation = (
            self.get_work_factor_search_given_nchildren_per_level(nchildren_per_level_simulation))


        work_factor_of_search_analytic = (
            self.get_work_factor_search_given_nchildren_per_level(nchildren_per_level_analytic))

        work_factor_permutations_analytic = (
            self.estimate_work_factor_permutations(attack_param_w, attack_param_la))

        work_factor_analytic = work_factor_of_search_analytic*work_factor_permutations_analytic
        work_factor_with_search_simulation = work_factor_of_search_simulation*work_factor_permutations_analytic

        data = {
            'av_size_of_LWSK': av_size_of_LWSK,
            'nchildren_per_level_simulation': nchildren_per_level_simulation,
            'nchildren_per_level_analytic': nchildren_per_level_analytic,
            'fraction_of_keys': fraction_of_keys,
            'work_factor_of_search_simulation': work_factor_of_search_simulation,
            'work_factor_of_search_analytic': work_factor_of_search_analytic,
            'work_factor_permutations_analytic': work_factor_permutations_analytic,
            'work_factor_analytic': work_factor_analytic,
            'work_factor_with_search_simulation': work_factor_with_search_simulation,
        }

        return data


def sparse_binary_matrix_to_file(filename, M):
    f = open(filename, 'w')

    f.write('%d %d\n' % (M.nrows(), M.ncols()))

    count = 0
    for i, row in enumerate(M):
        count += len(row.support())

    f.write('%d\n' % (count))

    for i, row in enumerate(M):
        for j, v in enumerate(row):
            if v:
                f.write('%d %d\n' % (i, j))

    f.close()


def secret_subset_to_file(filename, permutation):
    f = open(filename, 'w')

    f.write('%d\n' % len(permutation))

    for t in permutation:
        f.write('%d\n' % t)

    f.close()

def write_param_la(filename, la):
    f = open(filename, 'w')

    f.write('%d\n' % la)

    f.close()

def matrix_to_jcf_format(filename, M):
    f = open(filename, 'w')

    f.write('%d %d 2\n' % (M.nrows(), M.ncols()))

    count = 0
    for i, row in enumerate(M):
        count += len(row.support())

    f.write('%d\n' % (count))

    for i, row in enumerate(M):
        for j, v in enumerate(row.support()):
            if j == 0:
                f.write('-')
            f.write('%d\n' % (v + 1))

    f.close()

def get_A_low_weight_head(A, LWSA):

    head = list(LWSA)
    tail = list(A[len(head):])

    return matrix(GF(2), head + tail)

def find_secret_subset(Vsec, Vpub, LWSA, LWSK):
    pi = get_column_permutation_from_B_to_A(Vsec, Vpub)

    list_LWSK = list(LWSK)

    def perm_pi(v):
        print(v)
        return vector([v[pi[i]] for i in range(len(v))])

    secret_subset = [list_LWSK.index(perm_pi(v)) for v in LWSA]

    return secret_subset

def get_LWSA_linearly_independent(LWSA):
    LA = matrix(LWSA)

    col_echelon = LA.T.echelon_form()
    li_indexes = [c.support()[0] for c in col_echelon[:LA.rank()]]

    return matrix([LWSA[i] for i in li_indexes])

def gen_challenge(m, n, l, attack_param_w, attack_param_la, directory):

    if os.path.exists(directory):
        raise ValueError(f'Path {directory} already exists.')
    os.mkdir(directory)

    bpkp = BinaryPKP(m, n, l)
    A, Vsec, Vpub = bpkp.keygen_weak(attack_param_w, attack_param_la)
    attack = Attack(A, Vpub, attack_param_w, attack_param_la)

    LWSA_all, LWSK = attack.phase1_find_low_weight_keywords()
    LWSA = get_LWSA_linearly_independent(LWSA_all)[:attack_param_la]

    secret_subset_indexes = find_secret_subset(Vsec, Vpub, LWSA[:attack_param_la], LWSK)

    A_low_weight_head = get_A_low_weight_head(A, LWSA[:attack_param_la])

    matrix_to_jcf_format(os.path.join(directory, 'A_low_weight_head.jcf'), A_low_weight_head)
    matrix_to_jcf_format(os.path.join(directory, 'Vpub.jcf'), Vpub)
    sparse_binary_matrix_to_file(os.path.join(directory, 'LWSA.sparse_matrix'), matrix(LWSA))
    sparse_binary_matrix_to_file(os.path.join(directory, 'LWSK.sparse_matrix'), matrix(LWSK))
    secret_subset_to_file(os.path.join(directory, 'secret_subset.indexes'), secret_subset_indexes)
    write_param_la(os.path.join(directory, 'la.param'), attack_param_la)


    return A, Vsec, Vpub, LWSA, LWSK

def test_low_weight_sets_computation():
    bpkp = BinaryPKP(15, 38, 10)
    A, Vsec, Vpub = bpkp.keygen_weak(7, 8)
    attack = Attack(A, Vpub, 7, 8)
    LWSA, LWSK = attack.phase1_find_low_weight_keywords()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['generate_challenge', 'estimate_complexity'])
    parser.add_argument('m', type=int)
    parser.add_argument('n', type=int)
    parser.add_argument('l', type=int)
    parser.add_argument('attack_param_w', type=int)
    parser.add_argument('attack_param_la', type=int)
    parser.add_argument('--directory', help='Path to the directory where a challenge will be generated.')

    args = parser.parse_args()

    if args.action == 'generate_challenge':
        gen_challenge(args.m,
                      args.n,
                      args.l,
                      args.attack_param_w,
                      args.attack_param_la,
                      args.directory)

    elif args.action == 'estimate_complexity':
        estim = Estimates(args.m, args.n, args.l)
        estimates = estim.get_all_estimates_for_attack_parameters(args.attack_param_w,
                                                                  args.attack_param_la)
        pprint.pprint(estimates)

