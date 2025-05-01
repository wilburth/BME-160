# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A

'''
triadRewrite is a re-written version of the class Triad along with the same main() function from coordinateMathSoln.py

coordinateMathSoln asks for and collects three sets of coordinates using input() and will output the N-C and N-Ca bond lengths and the C-N-Ca bond angle.
Coordinates must be written on a single line.
    example input: C = (39.447, 94.657, 11.824) N = (39.292, 95.716, 11.027) Ca = (39.462, 97.101, 11.465)
'''

import math

class Triad:
    def __init__(self,p,q,r) :
        """ Construct a Triad.

        Example object construction:
            p1 = Triad( p=(1.,0.,0.), q=(0.,0.,0.), r=(0.,0.,0.) ).
        """
        self.p = p
        self.q = q
        self.r = r

        # self.pq1 = (p[0] - q[0])**2 + (p[1] - q[1])**2 + (p[2] - q[2])**2
        # self.pr1 = (p[0] - r[0])**2 + (p[1] - r[1])**2 + (p[2] - r[2])**2
        # self.qr1 = (q[0] - r[0])**2 + (q[1] - r[1])**2 + (q[2] - r[2])**2

        # instead of using zip() to iterate, do it manually --> [0], [1]. and [2] for (x,y,z)
        # bond length calculations without the sqrt
        # i.e. (ai - bi)^2 for i = x,y,z and a,b are the points
        self.pq2 = sum((p[i] - q[i])**2 for i in range(0, 3)) # from 0 up to but not including 3
        self.pr2 = sum((p[i] - r[i])**2 for i in range(0, 3))
        self.qr2 = sum((q[i] - r[i])**2 for i in range(0, 3))

        # self.angleP1 = (q[0] - p[0])*(r[0] - p[0]) + (q[1] - p[1])*(r[1] - p[1]) + (q[2] - p[2])*(r[2] - p[2])
        # self.angleQ1 = (p[0] - q[0])*(r[0] - q[0]) + (p[1] - q[1])*(r[1] - q[1]) + (p[2] - q[2])*(r[2] - q[2])
        # self.angleR1 = (p[0] - r[0])*(q[0] - r[0]) + (p[1] - r[1])*(q[1] - r[1]) + (p[2] - r[2])*(q[2] - r[2])

        # instead of using zip() to iterate, do it manually --> [0], [1], and [2] for (x,y,z)
        # bond angle calculations without the inverse cosign, without the division
        # i.e (ai - bi)*(ci - bi) for i = x,y,z and a,b,c are the points
        self.angleP2 = sum((q[i] - p[i])*(r[i] - p[i]) for i in range(0, 3)) # from 0 up to but not including 3
        self.angleQ2 = sum((p[i] - q[i])*(r[i] - q[i]) for i in range(0, 3))
        self.angleR2 = sum((p[i] - r[i])*(q[i] - r[i]) for i in range(0, 3))

    def dPQ (self):
        """ Provides the distance between point p and point q """
        # finish the bond length equation by adding the sqrt
        pq = math.sqrt(self.pq2)
        return pq

    def dPR (self):
        """ Provides the distance between point p and point r """
        pr = math.sqrt(self.pr2)
        return pr

    def dQR (self):
        """ Provides the distance between point q and point r """
        qr = math.sqrt(self.qr2)
        return qr

    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        # finish bond angle equation by adding inverse cosign, and dividing summation by
        angleP = math.acos(self.angleP2 / (math.sqrt(self.pq2) * math.sqrt(self.pr2)))
        return angleP

    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        angleQ = math.acos(self.angleQ2 / (math.sqrt(self.qr2) * math.sqrt(self.pq2)))
        return angleQ

    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        angleR = math.acos(self.angleR2 / (math.sqrt(self.pr2) * math.sqrt(self.qr2)))
        return angleR



def main():
    ''' main() will ask for and collect coordinates, then print out the N-C bond length, the N-Ca bond angle, and the C-N-Ca bond angle.'''

    triadInput = input('Enter Coordinates: ') # get triad sequence // example input: C = (39.447, 94.657, 11.824) N = (39.292, 95.716, 11.027) Ca = (39.462, 97.101, 11.465)
    triadCleaned = triadInput.replace('(', ',').replace(')', ',').split(',') #replace '(' and ')' with ',' then split at each ',' --> gives us our list!
    C = float(triadCleaned[1]), float(triadCleaned[2]), float(triadCleaned[3]) # sets 'C' equal to the stored data at points 1-3 in float format
    N = float(triadCleaned[5]), float(triadCleaned[6]), float(triadCleaned[7]) # sets 'N' equal to the stored data at points 5-7 in float format
    Ca = float(triadCleaned[9]), float(triadCleaned[10]), float(triadCleaned[11]) # sets 'Ca' equal to the stored data at points 9-11 in float format
    triadFinal = Triad(C, N, Ca)

    print('N-C bond length = {0:0.2f}'.format(triadFinal.dPQ())) # use method dPQ() to return the N-C bond length (to second demical)
    print('N-Ca bond length = {0:0.2f}'.format(triadFinal.dQR())) # use method dQR() to return the N-Ca bond length (to second decimal)
    print('C-N-Ca bond angle = {0:0.1f}'.format(triadFinal.angleQ() * (180.0 / math.pi))) # use method angleQ() to return C-N-Ca bond angle (to first decimal), along with changing the radian output to degrees

main()