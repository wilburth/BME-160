import sys

class FastAreader :
    '''
    Define the objects that are to be read by FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
    print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor then saves attribute fname '''
        self.fname = fname

    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta (self):
        ''' Read one entire FastA record. Then return the sequence header/sequence'''

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
#
# this is fully functional, but i would rlly appreciate if you tried to integrate this into your own program rather than copy/paste :)
# also i tried to describe why i did each thing the way i did so i hope that can also help you debug ur program too
#
# -----------GOALS---------------#
# methods to create power sets, find the uniques
# method to clean the sequence --> alternatively, can clean every sequence by hand but just easier to have a
# function do it for us
# method for printing output --> only use main to call methods from class FindUnique


class uniqueTRNA:
# try not to use global variables --> we can put everything inside of __init__ and instantiate it
    def __init__(self):
        # create power set list, a uniques list, and to make life easy, setup a dictionary which is
        # filled with the header + sequence from the FASTA file
        self.powerSets = []
        self.uniques = []
        self.headerSequenceDictionary = {}

        thisReader = FastAreader('someFileName.fa') # instantiation of fastareader with input file

        # to fill the headerSequenceDictionary, we can use a FOR LOOP along with a starting position of 0 to keep track
        position = 0
        for header, sequence in thisReader.readFasta(): # for each header + sequence that is read by our FastaReader
            cleanedSequence = self.clean(sequence) # call our .clean() method to clean the sequence that is provided from our FastaReader
            self.headerSequenceDictionary[position] = [header, cleanedSequence] # set each position in our headerSequenceDictionary to have a header as
            self.powerSets.append(self.powSet(cleanedSequence))
            position += 1


    def powSet(self, seq):
        # to find all of the unique sets, we first have to find all of the power sets that exist --> once we do this,
        # we can move onto the unique method which will sort through the power sets in order to find the true uniques that occur in the tRNA sequence
        powSet = set()
        for n in range(len(seq)): # for n in range of the length of the entire sequence
            temp = len(seq) # variable temp gets set to the length of the sequence as well --> we use this to keep track of where we are

            while temp > n: # while length2 is greater than POSITION n
                powSet.add(seq[n : temp]) # ADD a set to powSet which will be our cleanedSequence that has been CUT from POSITION n up to (but not uncluding) variable temp
                temp -= 1 # need to count down our temp variable by 1 during each loop in the while loop
                # this while loop exits when temp <= n

        return powSet

    def unique(self):
        # now that powSet has returned all of the existing powerSets within the tRna sequence, we can scan through and find the uniques
        # by creating duplicates of the powerSet list
        for sets in self.powerSets: # for all the existing sets that are in our powerSets list
            temp = set()
            copiedList = self.powerSets.copy() # create variable 'copiedList' which is a copy of the list of existing powerSets
            copiedSet = sets.copy() # create variable 'copiedSet' which is a copy of a SET in powerSets list
            copiedList.remove(copiedSet) # REMOVE the copiedSet (which is a single set) from the list of exisitng powerSets (this occurs ONLY if it exists)

            for powSet in copiedList: # for each power set in copiedList
                temp = temp.union(powSet) # use object 'union', call function union() with parameter powSet --> https://www.w3schools.com/python/ref_set_union.asp

            copiedSet.difference_update(temp) # compute difference on copySet with parameter as object 'union' --> https://www.w3schools.com/python/ref_set_difference_update.asp

            copiedSet2 = copiedSet.copy() # create variable copiedset2 which is a copy of a copy of the sets in powSetList with the difference_update method applied to it

            for tRna in copiedSet: # for each tRna that exists in copiedSet
                unique = copiedSet.copy() # create variable unique and set it to a copy of the copiedSet (which has had the difference update applied to it)
                unique.remove(tRna) # remove that tRna from the uniques

                for tRna2 in unique: # for all tRna that exist in unique
                    if tRna in tRna2 and len(tRna) < len(tRna2): # and if tRna1 is in tRna2 along with tRna #2 being longer
                        copiedSet2.discard(tRna2) # remove it from our copySet 2, and keep the essentials--> https://www.w3schools.com/python/ref_set_discard.asp

            self.uniques.append(copiedSet2)

    def clean(self, seq):
        # the clean method has 1 argument, seq, which it will clean by replacing all characters that do not belong, along with returning 'clean' which is the cleaned sequence
        clean = seq.replace('_', '').replace('-','').replace('.', '') # "They are not part of the sequence and must be removed. [-\_\.] are alignment characters."
        return clean


    def out(self):
        # print fasta input file in alphabetical order
        sortedDict = dict(sorted(self.headerSequenceDictionary.items(), key=lambda i:i[1]))
        for index in sortedDict:
            print(sortedDict[index][0], '\n', sortedDict[index][1])

            for pos in range(len(sortedDict[index][1])):
                for string in self.uniques[index]:
                    if string == sortedDict[index][1][pos: pos + len(string)]:
                        print(('.' * pos) + string)
def main():
    tRNA = uniqueTRNA()
    tRNA.unique()
    tRNA.out()

if __name__ == "__main__":
    main()

