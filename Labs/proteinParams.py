# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A

'''
proteinParams.py will provide the amino acid count, amino acid composition, theoretical pI, molar and mass extinction
coefficients, and molecular weight for a given protein sequence.

Example input:
    VLSPADKTNVKAAW

Example output:
    Number of Amino Acids: 14
    Molecular Weight: 1499.7
    molar Extinction coefficient: 5500.00
    mass Extinction coefficient: 3.67
    Theoretical pI: 9.88
    Amino acid composition:
        A = 21.43%
        C = 0.00%
        D = 7.14%
        E = 0.00%
        F = 0.00%
        G = 0.00%
        H = 0.00%
        I = 0.00%
        K = 14.29%
        L = 7.14%
        M = 0.00%
        N = 7.14%
        P = 7.14%
        Q = 0.00%
        R = 0.00%
        S = 7.14%
        T = 7.14%
        V = 14.29%
        W = 7.14%
        Y = 0.00%

'''
class ProteinParam :
    '''
    class ProteinParam provides the amino acid count, amino acid composition, theoretical pI, molar and mass extinction
    coefficients, and molecular weight for a given protein sequence.
    '''
    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like

    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        '''
        init() will parse through the user input, find the characters that match the keys from dictionary aa2mw,
        and will create a list from those characters. It will use this list to create the keys for a dictionary,
        and store the values for each key as the number of times the key occurs within the list.
        '''
        # A, C, D, E, F, G, H, I, L, K, M, N, P, Q, R, S, T, V, Y, W
        # set prot equal to whatever data is in "protein" (user input), make it of type string (str), make it uppercase
        prot = str(protein).upper()
        protList = []                           # creat empty list to be filled

        for char in prot:                       # for characters in prot (user input string.upper)
            if char in self.aa2mw.keys():       # if the character MATCHES a key in aa2mw
                protList.append(char)           # APPEND said character to the empty list (add it to the back)
        # print(protList)

        self.protClean = ''.join(protList).upper() # condense protList and make it uppercase (separator is '')

        self.aaTotal = 0        # initialize aaTotal at int(0)
        self.comp = {}          # initialize composition dictionary as empty
        # for loop counts and assigns value to keys based on how many times they occur within the cleaned protein string
        for key in self.aa2mw.keys():
            self.comp[key] = self.protClean.count(key) # set the key's value equal to the amount of times the key occurs within the protClean string (which is the cleaned user input string)
            # print(self.comp)

    def aaCount (self):
        '''
        aaCount() will parse through all keys in the dictionary created in init(), look for the values
        attached to each key, summate these values, and fill self.aaTotal (originally set at 0) with this summation.
        '''
        for key in self.comp.keys():             # for aa in composition keys
            self.aaTotal += self.comp[key]      # self.aaTotal is originally set at int(0), we add value attached to each key (this iterates through the ENTIRE list, and is technically adding the 0's)
        return self.aaTotal

    def pI (self):
        '''
        pI() solves for the pH that has a charge closest to 0, within 2 decimal places (i.e. 0.01).
        It does this through calling _charge_() on the specified pH, which will return the net charge
        of a given protein at the specified pH.
        '''
        maxpH = 14.0
        minpH = 0.0
        while (maxpH - minpH) > 0.01: # being while loop with max pH = 14.0 and min pH = 0.0 --> these variables change as the while loop iterates, the only time it will stop is when (maxpH - minpH) is LESS THAN 0.01
            midpoint = (maxpH + minpH)/2         # find the midpoint (we will either be setting the maxpH or minpH equal to this value)
            charge = self._charge_(midpoint)     # call method charge on the midpoint (first round is like calling charge on 7)
            if charge < 0:                  # if the calculated charge is less than 0 (negative)
                maxpH = midpoint                 # set the high point equal to the midpoint
            else:                           # if the calculated charge is greater than 0 (positive)
                minpH = midpoint                 # set the low point equal to the midpoint
        return midpoint      # this only returns mid when the while loop breaks (i.e. when high-low is less than 0.01, a 2 decimal place approximation


    def aaComposition (self) :
        '''
        aaComposition() returns the dictionary made in init(), which is created by parsing through the user input
        for characters that match the keys from aa2mw and creating a cleaned list, then using that cleaned list to
        fill an empty dictionary with keys that have matching values equal to the amount of times they occur in the
        user input sequence.
        '''
        return self.comp  # this is the dictionary that was made in _init_


    def _charge_ (self, pH):
        '''
        _charge_() follows the equations given in Lab03.ipynb, which are used to calculate the net charge
        of the user input sequence. It takes in a specified pH from the pI() method and calculates the net charge
        of the given protein at said specified pH value.
        '''
        summN = ((10**self.aaNterm) / (10**self.aaNterm + 10**pH))  # factor out the aaNterm amount from summation
        for key in self.aa2chargePos.keys():                        # iterate through all keys in aa2chargePos
            summN += self.comp[key] * ((10**self.aa2chargePos[key]) / (10**self.aa2chargePos[key] + 10**pH)) # self.comp[key] returns value attached to key (amount of times the key occurs)

        summC = ((10**pH) / (10**self.aaCterm + 10**pH))        # factor out the aaCterm amount from summation
        for key in self.aa2chargeNeg.keys():                    # iterate through all keys in aa2chargeNeg
            summC += self.comp[key] * ((10**pH) / (10**self.aa2chargeNeg[key] + 10**pH))

        return summN - summC # return solved summation of N terminus minus solved summation of C terminus

    def molarExtinction (self):
        '''
        molarExtinction() calculates the molar extinction coefficient of light absorbance at 280nm using Y, W, C content
        '''
        meCount = 0
        for key in self.aa2abs280:          # iterate over all keys in self.aa2abs280 (Y,W,C)
            meCount += self.aa2abs280[key] * self.comp[key]  # get value from key in aa2abs280, multiply value by value from key in comp
                                                            # (each value in comp is equal to the number of times the given key occured in the original input string)
        return meCount

    def massExtinction (self):
        '''
        massExtinction() calculates the mass extinction coefficient using the molecular weight and molar extinction coefficient
        '''
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        '''
        molecularWeight() returns 0.0 if the user input is empty, otherwise returns the molecular weight of the protein
        '''
        if self.aaTotal == 0:             # if aaTotal (the number of aa in input string) is false / empty, then the MW must be 0
            return 0.0
        else:
            currTotal = 0
            for key in self.protClean:                    # iterates through all keys
                currTotal += self.aa2mw[key] - self.mwH2O # value is equal to molecular weight of key minus the mw of H2O, summate
            return self.mwH2O + currTotal                 # return the mw of H2O along with the summation from prev line


# Please do not modify any of the following.  This will produce a standard output that can be parsed

import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        myAAcomposition = myParamMaker.aaComposition()
        keys = list(myAAcomposition.keys())
        keys.sort()
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present
        for key in keys :
            print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber))

        inString = input('protein sequence?')

if __name__ == "__main__":
    main()

