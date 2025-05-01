# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A

"""
sequenceAnalysis Module docstring:

sequenceAnalysis is a module with classes OrfFiner, NucParams, FastAreader, and ProteinParams, each of which can be
imported into a python script. Refer to each docstring for their respective purposes.
"""

# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A

import numpy as np
class Hydrophobicity:
    """
    class docstring:

    """
    rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'   # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    hydropathyTable = {
        # hydropathy scores
        # positive = hydrophobic
        'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
        # negative = hydrophilic
        'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
        'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9,
        'R': -4.5,
        # nothing value (codes for no AA)
        '-': '-STOP-'
    }

    def __init__(self, seq, head, frame, start, stop, length):
        """
        __init__() docstring:

        """
        self.seq = seq
        self.head = head
        self.frame = str(frame)
        self.start = str(start)
        self.stop = str(stop)
        self.length = str(length)

        self.codonList = [] # create our codons from input seq
        self.codonToAA = [] # match codons to amino acids
        self.aaToValue = [] # match amino acids to hydropathy value

        self.codonList = [self.seq[i:i+3] for i in range(0, len(self.seq), 3)]

        for item in self.codonList:
            if item in self.rnaCodonTable.keys():
                matchingAA = self.rnaCodonTable[item]
                self.codonToAA.append(matchingAA)

        for item in self.codonToAA:
            if item in self.hydropathyTable.keys():
                matchingValue = self.hydropathyTable[item]
                self.aaToValue.append(matchingValue)

    #https://stackoverflow.com/questions/14313510/how-to-calculate-rolling-moving-average-using-numpy-scipy
    #def moving_average(self, w):
        #return np.convolve(self, np.ones(w), 'valid') / w


    def out(self, outFile):

        outFile.write(self.head + '\n')
        outFile.write('Frame: ' + self.frame + '\n')
        outFile.write('Starting Pos: ' + self.start + '\n')
        outFile.write('Ending Pos: ' + self.stop + '\n')
        outFile.write('Length: ' + self.length + '\n')
        outFile.write('Gene Sequence: ' + self.seq + '\n')

        for i in range(0, len(self.codonList)):
            output1 = '{0} . . . {1} . . . {2} . . . {3}\n'.format(i ,self.codonList[i], self.codonToAA[i], self.aaToValue[i])
            outFile.write(output1)

        outFile.write('\n')


# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A

class OrfFinder:
    """
    Class OrfFinder has methods __init__(), orf(), compliment(), processingData(), and output(), which work together to
    parse through a genome and find Open Reading Frames (ORFs).
    """

    def __init__(self, seq, header, start, stop, minGene, longest):
        """
        method __init__() will declare the default start, stop, minimum gene length, and longest open reading frame,
        unless otherwise stated on the commandLine by the user.
        """
        self.seq = seq
        self.header = header

        if start:
            self.start = start
        else:
            self.start = ['ATG']

        if stop:
            self.stop = stop
        else:
            self.stop = ['TAG', 'TGA', 'TAA']

        if minGene:
            self.min = minGene
        else:
            self.min = 100

        if longest:
            self.longest = True
        else:
            self.longest = False


        # self.start = ['ATG'] # default start
        # self.stop = ['TAG','TGA','TAA'] # default stop
        # self.min = 100 # default min gene length
        # self.longest = True # default for longest orf gene

        self.orfData = []



    def orf(self, compliment=False):
        """
        method orf() will parse through and find the codons in the open reading frames.

        """
        #                                123
        # 6 total reading frames --> 5'- ATG AAA AAA -3'
        #                            3'- AAA AAA GTA -5'    (complementary strand)
        #                                      -(321)
        #
        # focusing just on frames 1,2,3 --> can reverse / flip them and re-run orf --> need def reverse() / compliment()
        # have to keep track of start, stop and orientation
        # (start position, stop position, frame orientation)

        # if complimentary strand is true:
        #   reverse the sequence
        #   set orientation to -1
        # else (complimentary strand is false):
        #   set orientation to 1
        #
        # for frame in range3 (internal 0,1,2 translates to external 1,2,3)
        #   set frame orientation
        #   for i in range(frames, length of sequence, split of 3)
        #       create codons --> .join()
        #       if codon in start
        #           append to startList
        #       if codon in stop
        #           record stop pos
        #           for start in startList
        #               call class(start, stop, frame orientation) to work with items

        startLocations = [0] # first element of list used to handle dangling start case --> +1 1...393 393 doesn't catch without this

        # seq = self.seq
        # orientation = 1 # default orientation until figure out how to do the -1,-2,-3

        if compliment:
            # [::-1] to reverse the sequence
            # maketrans() to swap characters
            # use translate() --> https://python-reference.readthedocs.io/en/latest/docs/str/translate.html --> example 2
            seq = self.seq.translate(str.maketrans('ATGCUN', 'TACGAN'))[::-1]
            orientation = -1
        else:
            seq = self.seq
            orientation = 1

        for frame in range(3): # going through (external) frames(1,2,3) using (internal) position(0,1,2)
            frameOrientation = (frame+1) * orientation # need frame+1 since pos(0,1,2) maps to frame(1,2,3)

            for i in range(frame, len(seq), 3): # for i in range(starting at frame, up to length of sequence, step of 3)
                codon = ''.join(seq[i:i+3]) # creates our codons

                if codon in self.start: # if codon matches start parameter --> append its LOCATION to a list
                    startLocations.append(i)  # location appended in startLocations list

                if codon in self.stop: # if codon matches stop parameter --> get its stop position (need +2), send it to another method to process it
                    stopPosition = i + 2

                    if self.longest and len(startLocations) > 0: # this if loop runs only the largest ORF in the gene
                        startLocations = [startLocations[0]] # assume largest orf is the first element in startLocations (sorted this way)

                    for start in startLocations: # for each start position in startLocations list
                        self.processingData(start, stopPosition, frameOrientation)    # call processingData with parameters of start, stop, and orientation
                    startLocations = [] # clear the list

                if i is len(seq):   # for i is the length of the sequence --> dangling stop case
                    stopPosition = i    # set stop position to last position

                    for start in startLocations:    # for all starts in startLocations list
                        self.processingData(start, stopPosition, frameOrientation)  # call processingData with parameters start, stop --> has been overwritten with i as the position instead of i+2, and orientation
                    startLocations = [] # clear the list


    def compliment(self):
        """
        helper function for orf(). Will rerun orf() method with compliment as True value, which sets orientation to be
        reverse (starting at -1).

        """
        self.orf(compliment=True)


    def processingData(self, start, stop, orientation):
        """
        helper function for orf(). Will process incoming start, stop, and frame orientation data from orf(),
        and append it to the main orfData list.
        """
        # stop should always be greater than start --> need to make sure this works when flipping sequence for complementary seq
        # everything is 1 off ??? --> +1 fix


        sequenceLength = stop - start + 1

        if orientation < 0: # handles frame cases for reverse strand
            temp = stop
            stop = len(self.seq) - start - 1
            start = len(self.seq) - temp - 1

        if sequenceLength >= self.min:  # make sure sequence length is greater than parameter self.min (default=100)
            item = [orientation, start + 1, stop + 1, sequenceLength] # creating a list as "item"
            self.orfData.append(item)  # append list to list orfData --> array 2D


    def output(self, outFile):
        '''
        outputFile sorts dataList by order of largest length and then by leftmost gene location in the case of a tie
        and then writes the sorted list to an output file.
        '''
        outFile.write('{0}\n'.format(self.header))
        for item in sorted(self.orfData, key=lambda i: (i[3], i[1]), reverse=True):   # sorts orf data list by decreasing seq length, then by left position of the gene
            output = '{:+d} {:>5d}..{:>5d} {:>5d}\n'.format(item[0], item[1], item[2], item[3]) # formatted print (frame, leftgene, rightgene, eqlength)
            outFile.write(output)




class NucParams:
    """
    class NucParams docstring:

    class NucParams builds 3 dictionaries (aaComp, nucComp, codonComp) in the __init__ method. The addSequence method
    parses through an input sequence and fills these dictionaries with values that are based on how many times they
    occur within the input sequence (denoted as parameter 'inSeq').
    """
    rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString = ''):
        """
        __init__() docstring:

        __init__() will initialize all necessary dictionaries to be used by genomeAnalyzer.py.
        The aaComp dictionary will have keys that are amino acids, and the corresponding values will be how many times
        the amino acid is coded within the sequence.
        The nucComp dictionary will have keys that are A,T,G,C,U,N, and the corresponding values will be how many times
        the nucleotide appears within the sequence.
        The codonComp dictionary will have keys that are codons, and the corresponding values will be how many times t
        the codon appears within the sequence.
        """

        # dictionary used to count amino acids
        # has the amino acid as the key with 0 as the initial value
        self.aaComp = {}
        for value in self.rnaCodonTable.values():
            self.aaComp[value] = 0

        # dictionary used to count individual nucleotides, only counting ACGTUN (no N in testGenome.fa ?)
        # has the nucleotide as the key with 0 as the initial value
        allowedChar = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'U': 5, 'N': 6}
        self.nucComp = {}
        for key in allowedChar.keys():
            self.nucComp[key] = 0

        # dictionary used to count codons (triplets)
        # has the codon as the key with 0 as the initial value
        self.codonComp = {}
        for key in self.rnaCodonTable.keys():
            self.codonComp[key] = 0

        self.addSequence(inString)

    def addSequence (self, inSeq):
        """
        addSequence() docstring:

        addSequence() will parse through inSeq and fill the dictionaries made in _init_ (aaComp, nucComp, codonComp).

        """

        # for character in inSeq (user input string, passed from _init_)
        for char in inSeq:
            if char in self.nucComp:             # if the character matches a nucleotide
                self.nucComp[char] += 1          # add 1 count to the characters value

        rnaSeq = inSeq.replace('T', 'U') # change the DNA strand into RNA strand (necessary)

        # following for / if loop is based on this idea from stackoverflow
        # self.y = [self.x[i:i+3] for i in range(0, len(self.x), 3)]
        # https://stackoverflow.com/questions/9475241/split-string-every-nth-character
        #
        # only have to parse through self.codonComp.keys to fill out the codonComp and aaComp tables, because we will
        # never be adding a value to an aaComp key unless that key is also a value in codonComp.values()

        codonList = [rnaSeq[i:i + 3] for i in range(0, len(rnaSeq), 3)] # working as intended

        for item in codonList:                  # for item (position) in codonList
            if item in self.codonComp.keys():   # if the item is also a key in dictionary codonComp
                self.codonComp[item] += 1       # add 1 count to the value of the key
                matchingAmino = self.rnaCodonTable[item]        # can also count aa's in this loop --> set aa to be what codon codes for
                self.aaComp[matchingAmino] += 1                 # add 1 count to the value of the aa in dictionary aaComp


    def aaComposition(self):
        """
        aaComposition() docstring:

        Returns returns the filled amino acid dictionary (aaComp) to be used in genomeAnalyzer.py
        """
        return self.aaComp

    def nucComposition(self):
        """
        nucComposition() docstring:

        Returns the filled nucleotide dictionary (nucComp) to be used in genomeAnalyzer.py
        """
        return self.nucComp

    def codonComposition(self):
        """
        codonComposition() docstring:

        Returns the filled codon dictionary to be used in genomeAnalyzer.py
        """
        return self.codonComp

    def nucCount(self):
        """
        nucCount() docstring:

        Returns the sum of all values in the nucleotide dictionary (total # of nucleotides)
        to be used in genomeAnalyzer.py
        """
        return sum(self.nucComp.values())


# for integration with findUnique.py, class FastAreader was updated by Dr. B in the Lab06.ipynb file.
# the change is an addition of lines, which read 'if not line:', followed by 'return header, sequence'
import sys
class FastAreader :

    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta (self):

        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                if not line: # we are at EOF
                    return header, sequence
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence


class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog -  ',
                                              epilog = 'Program epilog - ',
                                              add_help = True, #default is True
                                              prefix_chars = '-',
                                              usage = '%(prog)s <infile> <outfile> <outfile2> <option?s> || example usage: findORFs.py inFile.fa outFile.txt -mG=XXX -lG'
                                              )
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        self.parser.add_argument('outFile', action = 'store', help='output file name')

        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False,
                                 help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (100,200,300,500,1000), default=100,
                                 action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?',
                                 help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?',
                                 help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)




# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A


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




