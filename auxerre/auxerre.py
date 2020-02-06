"""
auxerre.py

Transport Properties from Molecular Dynamics Using Reciprocal-Space Correlation Analysis

_inteload: http://mdtraj.org/development/api/generated/mdtraj.iterload.html

"""

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import itertools


class Calculator(object):
    """
    A base class for all transport property calculators.

    """
    def process_trajectory(self, trajectory, cosines, sines):
        pass


class DiffusivityCalculator(Calculator):
    """
    Diffusion coefficient calculator.

    """
    def __init__(self, cross_correlations_only=False):
        self._cross_correlations_only = cross_correlations_only
        self._A = []
        self._B = []

    def digest_trajectory(self, trajectory, cosines, sines):
        self._A.append(np.sum(cosines, axis=2))
        self._B.append(np.sum(sines, axis=2))

    def correlate(self):
        A = np.concatenate(self._A, axis=1)
        B = np.concatenate(self._B, axis=1)
        for a, b in zip(A, B):
            # print(a, b)
            v = np.correlate(a, a, mode='same') + np.correlate(b, b, mode='same')
            plt.plot(v)
        # print(A.shape)
        plt.show()

class Analyzer(object):
    """
    Molecular Dynamics trajectory analyzer.

    Parameters
    ----------
        calculators : list(:class:`Calculator`)
            A list of :class:`Calculator` objects.
        box_lengths : list or np.array
            The lengths of the simulation box in nanometers.

    Keyword Args
    ------------
        cutoff : int, default=5
            The maximum length of the integer lattice vectors considered in the analysis.

    """
    def __init__(self, calculators, box_lengths, cutoff=5):
        if isinstance(calculators, Calculator):
            self._calculators = [calculators]
        else:
            self._calculators = calculators
        box_vectors = np.diagflat(np.array(box_lengths))

        self._q_vectors = self._generate_q_vectors(cutoff, box_vectors)
        self._transposed_q_vectors = np.transpose(self._generate_q_vectors(cutoff, box_vectors))

    def _generate_q_vectors(self, cutoff, box_vectors):
        """
        """
        M = 2*cutoff + 1    # Number of integers from -cutoff to +cutoff
        maxnvecs = M**3//2  # Cube of odd is odd -> division by 2 rounds down
        MM = M*M
        vectors = []
        squared_lengths = []
        squared_cutoff = cutoff*cutoff
        for index in range(maxnvecs):
            i = index + maxnvecs + 1
            n1 = i//MM
            j = i - n1*MM
            n2 = j//M
            n3 = j - n2*M
            vector = np.array([n1, n2, n3], dtype=np.int) - cutoff
            squared_length = np.dot(vector, vector)
            if squared_length <= squared_cutoff:
                vectors.append(vector)
                squared_lengths.append(squared_length)
        lattice_vectors = [vectors[i] for i in np.argsort(squared_lengths)]
        reciprocal_matrix = np.linalg.inv(box_vectors)
        q_vectors = 2.0*np.pi*np.array(lattice_vectors).dot(reciprocal_matrix)
        return q_vectors

    def add_trajectory(self, trajectory):
        # args = np.einsum('ij,klj->ikl', self._q_vectors, trajectory.xyz)
        print('digesting trajectory...')
        args = np.einsum('ijk,lk->lij', trajectory.xyz, self._q_vectors)
        cosines = np.cos(args)
        sines = np.sin(args)
        for calculator in self._calculators:
            calculator.digest_trajectory(trajectory, cosines, sines)

    def correlate(self):
        for calculator in self._calculators:
            calculator.correlate()
