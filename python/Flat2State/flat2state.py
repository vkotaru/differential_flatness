import numpy as np 
import scipy

def hat(vector):
    return np.array([[0., -vector[2], vector[1]], 
                        [vector[2], 0., -vector[0]],
                        [-vector[1], vector[0], 0.]])

def vee(matrix):
    return np.array([matrix[2,1], matrix[0,2], matrix[1,0]])


class Flat2State(object):
    def __init__(self):
        pass 

    @staticmethod
    def compute_moment(axQ,  daxQ, d2axQ, mQ=0.85, J=None):
        # variables
        if J is None:
            J = np.array([[0.005315307431627, 0.000005567447099,  0.000005445855427],
                            [0.000005567447099, 0.004949258422243, 0.000020951458431],
                            [0.000005445855427, 0.000020951458431, 0.009806225007686]])
        g = 9.80655
        e1 = np.array([1.0,0.0,0.0])
        e2 = np.array([0.0,1.0,0.0])
        e3 = np.array([0.0,0.0,1.0])   

        b1d = e1
        db1d = np.zeros((3))
        d2b1d = np.zeros((3)) 

        fb3 = np.dot(mQ, axQ+np.dot(g, e3))
        norm_fb3 = np.linalg.norm(fb3)
        f = norm_fb3
        b3 = np.divide(fb3, norm_fb3)
        b3_b1d = np.cross(b3, b1d)
        norm_b3_b1d = np.linalg.norm(b3_b1d)
        b1 = np.divide(-np.cross(b3, b3_b1d), norm_b3_b1d)
        b2 = np.cross(b3, b1)
        R = np.array(np.vstack((b1, b2, b3)))

        dfb3 = np.dot(mQ, daxQ)
        dnorm_fb3 = np.divide(fb3.dot(dfb3), norm_fb3)
        db3 = np.divide(np.multiply(dfb3, norm_fb3)-np.multiply(fb3, dnorm_fb3), norm_fb3**2.)
        db3_b1d = np.cross(db3, b1d)+np.cross(b3, db1d)
        dnorm_b3_b1d = np.divide(b3_b1d.dot(db3_b1d), norm_b3_b1d)
        db1 = np.divide(-np.cross(db3, b3_b1d)-np.cross(b3, db3_b1d)-np.multiply(b1, dnorm_b3_b1d), norm_b3_b1d)
        db2 = np.cross(db3, b1)+np.cross(b3, db1)
        dR = np.array(np.vstack((db1, db2, db3)))

        Omega = vee(np.dot(dR, R.conj().T))

        d2fb3 = np.dot(mQ, d2axQ)
        d2norm_fb3 = np.divide(dfb3.dot(dfb3)+fb3.dot(d2fb3)-np.dot(dnorm_fb3, dnorm_fb3), norm_fb3)
        d2b3 = np.divide(np.dot(np.dot(d2fb3, norm_fb3)+np.dot(dfb3, dnorm_fb3)-np.dot(dfb3, dnorm_fb3)-np.dot(fb3, d2norm_fb3), norm_fb3**2.)-np.dot(np.dot(np.dot(db3, norm_fb3**2.)*2., norm_fb3), dnorm_fb3), norm_fb3**4.)
        d2b3_b1d = np.cross(d2b3, b1d)+np.cross(db3, db1d)+np.cross(db3, db1d)+np.cross(b3, d2b1d)
        d2norm_b3_b1d = np.divide(np.dot(db3_b1d.dot(db3_b1d)+b3_b1d.dot(d2b3_b1d), norm_b3_b1d)-np.dot(b3_b1d.dot(db3_b1d), dnorm_b3_b1d), norm_b3_b1d**2.)
        d2b1 = np.divide(np.dot(-np.cross(d2b3, b3_b1d)-np.cross(db3, db3_b1d)-np.cross(db3, db3_b1d)-np.cross(b3, d2b3_b1d)-np.dot(db1, dnorm_b3_b1d)-np.dot(b1, d2norm_b3_b1d), norm_b3_b1d)-np.dot(np.dot(db1, norm_b3_b1d), dnorm_b3_b1d), norm_b3_b1d**2.)
        d2b2 = np.cross(d2b3, b1)+np.cross(db3, db1)+np.cross(db3, db1)+np.cross(b3, d2b1)
        d2R = np.array(np.vstack((d2b1, d2b2, d2b3)))
        dOmega = vee((np.dot(dR, dR.conj().T)+np.dot(d2R, R.conj().T)))
        #%vee( dR'*dR + R'*d2R, true ) ;
        M = np.dot(J, dOmega)+np.cross(Omega, np.dot(J, Omega))

        s = {}
        s['R'] = R
        s['Omega'] = Omega
        s['dOmega'] = dOmega
        s['M'] = M
        s['f'] = f

        return s

    @staticmethod
    def quadrotor(traj, mQ = .85, J = None):
        # variables
        if J is None:
            J = np.array([[0.005315307431627, 0.000005567447099,  0.000005445855427],
                            [0.000005567447099, 0.004949258422243, 0.000020951458431],
                            [0.000005445855427, 0.000020951458431, 0.009806225007686]])

        g = 9.80655
        e1 = np.array([1.0,0.0,0.0])
        e2 = np.array([0.0,1.0,0.0])
        e3 = np.array([0.0,0.0,1.0])    

        # flats
        xQ = traj['x']
        vxQ = traj['dx']
        axQ = traj['d2x']
        daxQ = traj['d3x']
        d2axQ = traj['d4x']

        s = Flat2State.compute_moment(axQ, daxQ, d2axQ, mQ, J)   
        s['xQ'] = xQ
        s['vQ'] = vxQ
        s['aQ'] = axQ

        return s
