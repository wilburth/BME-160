# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A

'''
coordinateMathSoln asks for and collects three sets of coordinates using input() and will output the N-C and N-Ca bond lengths and the C-N-Ca bond angle.
Coordinates must be written on a single line.
    example input: C = (39.447, 94.657, 11.824) N = (39.292, 95.716, 11.027) Ca = (39.462, 97.101, 11.465)

'''

import math
class Triad :
    """
    Calculate angles and distances among a triad of points.

    Author: David Bernick
    Date: March 21, 2013
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()

    Required Modules: math
    initialized: 3 positional tuples representing Points in n-space
             p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,1,0) )
    attributes: p,q,r the 3 tuples representing points in N-space
    methods:  angleP(), angleR(), angleQ() angles measured in radians
          dPQ(), dPR(), dQR() distances in the same units of p,q,r

    """

    def __init__(self,p,q,r) :
        """ Construct a Triad.

        Example construction:
            p1 = Triad( p=(1.,0.,0.), q=(0.,0.,0.), r=(0.,0.,0.) ).
        """
        self.p = p
        self.q = q
        self.r = r



    # private helper methods
    def d2 (self,a,b) : # calculate squared distance of point a to b
        return float(sum((ia-ib)*(ia-ib)  for  ia,ib in zip (a,b)))

    def dot (self,a,b) : # dotProd of standard vectors a,b
        return float(sum(ia*ib for ia,ib in zip(a,b)))

    def ndot (self,a,b,c) : # dotProd of vec. a,c standardized to b
        return float(sum((ia-ib)*(ic-ib) for ia,ib,ic in zip (a,b,c)))

    # calculate lengths(distances) of segments PQ, PR and QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return math.sqrt(self.d2(self.p,self.q))

    def dPR (self):
        """ Provides the distance between point p and point r """
        return math.sqrt(self.d2(self.p,self.r))

    def dQR (self):
        """ Provides the distance between point q and point r """
        return math.sqrt(self.d2(self.q,self.r))

    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return math.acos(self.ndot(self.q,self.p,self.r) /   math.sqrt(self.d2(self.q,self.p)*self.d2(self.r,self.p)))

    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return math.acos(self.ndot(self.p,self.q,self.r) /  math.sqrt(self.d2(self.p,self.q)*self.d2(self.r,self.q)))

    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return math.acos(self.ndot(self.p,self.r,self.q) /  math.sqrt(self.d2(self.p,self.r)*self.d2(self.q,self.r)))

def main():
    ''' main() will ask for and collect coordinates, then print out the N-C bond length, the N-Ca bond angle, and the C-N-Ca bond angle.'''

    userInput = input('Enter Coordinates: ') # get triad sequence
    triadCleaned = userInput.replace('(', ',').replace(')', ',').split(',') #replace '(' and ')' with ',' then split at each ',' --> gives us our list!
    C = float(triadCleaned[1]), float(triadCleaned[2]), float(triadCleaned[3]) # set 'C' equal to the stored data at points 1-3 in float format
    N = float(triadCleaned[5]), float(triadCleaned[6]), float(triadCleaned[7]) # set 'N' equal to the stored data at points 5-7 in float format
    Ca = float(triadCleaned[9]), float(triadCleaned[10]), float(triadCleaned[11]) # set 'Ca' equal to the stored data at points 9-11 in float format
    triadFinal = Triad(C, N, Ca)

    print('N-C bond length = {0:0.2f}'.format(triadFinal.dPQ())) # use method dPQ() to return the N-C bond length (to second demical)
    print('N-Ca bond length = {0:0.2f}'.format(triadFinal.dQR())) # use method dQR() to return the N-Ca bond length (to second decimal)
    print('C-N-Ca bond angle = {0:0.1f}'.format(triadFinal.angleQ() * (180.0 / math.pi))) # use method angleQ() to return C-N-Ca bond angle (to first decimal), along with changing the radian output to degrees

main()