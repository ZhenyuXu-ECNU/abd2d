import numpy as np
from body import Body
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix


class OrthogonalEnergy:

    @staticmethod
    def val(body):
        a11 = body.q[2]
        a12 = body.q[3]
        a21 = body.q[4]
        a22 = body.q[5]
        a1 = np.array([a11, a12])  # row vec
        a2 = np.array([a21, a22])  # row vec
        coeff = body.volume * body.volume

        v = coeff * ((np.dot(a1, a1) - 1) ** 2 + (np.dot(a2, a2) - 1) ** 2)
        v += 2 * coeff * np.dot(a1, a2) ** 2
        return v

    @staticmethod
    def grad(body):
        a11 = body.q[2]
        a12 = body.q[3]
        a21 = body.q[4]
        a22 = body.q[5]
        a1 = np.array([a11, a12])  # row vec
        a2 = np.array([a21, a22])  # row vec
        coeff = body.volume * body.volume

        dvda1 = 2 * coeff * (2 * (np.dot(a1, a1) - 1) * a1 + 2 * np.dot(a2 * a2[:, None], a1))
        dvda2 = 2 * coeff * (2 * (np.dot(a2, a2) - 1) * a2 + 2 * np.dot(a1 * a1[:, None], a2))
        return np.array([0, 0, dvda1[0], dvda1[1], dvda2[0], dvda2[1]])

    @staticmethod
    def hess(body):
        a11 = body.q[2]
        a12 = body.q[3]
        a21 = body.q[4]
        a22 = body.q[5]
        a1 = np.array([a11, a12])  # row vec
        a2 = np.array([a21, a22])  # row vec
        coeff = body.volume * body.volume

        hvda1da1 = 2 * coeff * (4 * a1 * a1[:, None] + 2 * (np.dot(a1, a1) - 1) * np.identity(2) + 2 * (a2 * a2[:, None]))
        hvda2da2 = 2 * coeff * (4 * a2 * a2[:, None] + 2 * (np.dot(a2, a2) - 1) * np.identity(2) + 2 * (a1 * a1[:, None]))

        hvda1da2 = 2 * coeff * np.array([[4 * a11 * a21 + 2 * a12 * a22, 2 * a12 * a21],
                                         [2 * a11 * a22, 2 * a11 * a21 + 4 * a12 * a22]])

        return np.array([[0, 0,              0,              0,              0,              0],
                         [0, 0,              0,              0,              0,              0],
                         [0, 0, hvda1da1[0, 0], hvda1da1[0, 1], hvda1da2[0, 0], hvda1da2[0, 1]],
                         [0, 0, hvda1da1[1, 0], hvda1da1[1, 1], hvda1da2[1, 0], hvda1da2[1, 1]],
                         [0, 0, hvda1da2[0, 0], hvda1da2[1, 0], hvda2da2[0, 0], hvda2da2[0, 1]],
                         [0, 0, hvda1da2[0, 1], hvda1da2[1, 1], hvda2da2[1, 0], hvda2da2[1, 1]]])


class InertiaEnergy:

    def __init__(self, body):
        self.body = body

    @staticmethod
    def val(body):
        q       = body.q
        q_tilde = body.q_tilde
        M       = body.M
        return 0.5 * np.dot(q - q_tilde, np.dot(M, (q - q_tilde)))

    @staticmethod
    def grad(body):
        q       = body.q
        q_tilde = body.q_tilde
        M       = body.M
        return np.dot(M, q - q_tilde)

    @staticmethod
    def hess(body):
        M       = body.M
        return M


class TotalIpc:

    @staticmethod
    def val(body_list: list[Body]):
        res = 0

        for body in body_list:
            res += OrthogonalEnergy.val(body) + InertiaEnergy.val(body)

        return res

    @staticmethod
    def search_dir(body):
        H_body = csr_matrix(OrthogonalEnergy.hess(body) + InertiaEnergy.hess(body))
        grad = OrthogonalEnergy.grad(body) + InertiaEnergy.grad(body)
        return spsolve(H_body, -grad)



