# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A

"""
findORFs.py is a program intended to find the Open Reading Frames in a given genome, and return an organized
output file to the user which will contain the start, stop, length, and frame orientation of the open reading frame.

This program also implements the class CommandLine, which allows the user access to specific commands that can assist
in finding ORFs of a specific length, or those with only a specific start or stop codon.

This program outputs a file to the directory that the program is placed in -- the name will be what the user sets it to
be in the command line. Use -h or --help for more information on the usable commandLine arguments.
"""

from sequenceAnalysis import FastAreader
from sequenceAnalysis import OrfFinder

########################################################################
# CommandLine
########################################################################
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
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does',
                                              epilog = 'Program epilog - some other stuff you feel compelled to say',
                                              add_help = True, #default is True
                                              prefix_chars = '-',
                                              usage = '%(prog)s <infile> <outfile> <option?s> || example usage: findORFs.py inFile.fa outFile.txt -mG=XXX -lG'
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


########################################################################
# Main
# Here is the main program
#
#
########################################################################


def main(inCL=None):
    '''
    Find some genes.
    '''
    if inCL is None:
        myCommandLine = CommandLine() #['tass2.fa','tass2ORFdata-ATG-100.txt','--longestGene']
    else :
        myCommandLine = CommandLine(inCL)

    ###### replace the code between comments.
    #print (myCommandLine.args)

    longestGene = myCommandLine.args.longestGene # is True if only the longest Gene is desired
    start = myCommandLine.args.start # is a list of start codons
    minGene = myCommandLine.args.minGene # is the minimum Gene length to include
    stop = myCommandLine.args.stop # optional stop arguments

    orfReader = FastAreader(myCommandLine.args.inFile)  # runs FastAreader on command line infile
    with open(myCommandLine.args.outFile, 'w') as outFile: # opens command line outfile for writing
        for head, seq in orfReader.readFasta():

            orf = OrfFinder(seq, head, start, stop, minGene, longestGene)
            orf.orf()
            orf.compliment()
            orf.output(outFile)

    # myCommandLine.args.inFile has the input file name
    # myCommandLine.args.outFile has the output file name

    #for head, seq in orfReader.readFasta():
    #    orf = OrfFinder(seq, head, start, stop, minGene, longestGene)
    #    orf.orf()
    #    orf.output()




#######

if __name__ == "__main__":
    main()  # delete the list when you want to run with STDIN

