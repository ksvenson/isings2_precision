import numpy as np
from quat import Quat


class GrpElemO4:

    def __init__(self, ql, qr, inv):
        if (isinstance(ql, Quat)):
            self.ql = ql
        elif (isinstance(ql, np.ndarray) and ql.size == 4):
            self.ql = Quat(ql[0], ql[1], ql[2], ql[3])
        else:
            raise TypeError(f'Invalid quaternion parameter {ql}')

        if (isinstance(qr, Quat)):
            self.qr = qr
        elif (isinstance(qr, np.ndarray) and qr.size == 4):
            self.qr = Quat(qr[0], qr[1], qr[2], qr[3])
        else:
            raise TypeError(f'Invalid quaternion parameter {qr}')

        if (inv == 1 or inv == -1):
            self.inv = inv
        else:
            raise TypeError(f'Invalid inversion parameter {inv}')

    def __str__(self):
        return "{} {} {}".format(self.ql, self.qr, self.inv)

    def __mul__(self, rhs):
        if isinstance(rhs, GrpElemO4):
            # multiply two O(4) group elements
            if self.inv == -1:
                # for quaternions, (ab)^* = b^* a^* because they're non-abelian
                ql = self.ql * rhs.qr.Inverse()
                qr = rhs.ql.Inverse() * self.qr
            else:
                ql = self.ql * rhs.ql
                qr = rhs.qr * self.qr

            return GrpElemO4(ql, qr, self.inv * rhs.inv)

        elif isinstance(rhs, np.ndarray) and rhs.size == 4:
            # rotate a 4-vector by this group element
            q = Quat(rhs[3], rhs[0], rhs[1], rhs[2])
            if self.inv:
                q = q.Inverse()
            q = self.ql * q * self.qr
            return np.array([q.x, q.y, q.z, q.w])

        else:
            raise TypeError(f'Can\'t multiply GrpElemO4 by {rhs}')

    def Identity():
        return GrpElemO4(Quat(1.0, 0.0, 0.0, 0.0), Quat(1.0, 0.0, 0.0, 0.0), 1)
