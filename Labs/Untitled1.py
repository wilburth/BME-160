#!/usr/bin/env python
# coding: utf-8

# In[4]:


# Name: Andrea Ramos Coronado (aramosco)
# Group Members:None
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

class uniqueTRNA:
    '''
    The class collect the header and sequences from one FastA file then finds unique and essentials subsequences.
    The unique are found by creating new sets. New sets that only contain items that agrument 1 and agrument 2 don't have can be made by using the ^ function.
    Essentials are found by creating a list of the uniques and then iterating every agrument within uniques through the list. When a substring within uniques matches the an item within the list, then it's stored within a set.
    That set is later ^ from the unique sent to create a set only containing unique essentials
    Input: FastA file
    Ouput:
    Header
    sequence
    '.'for every none uniqueessential until it reaches an uniqueessential then it prints it out.
    '''
    TRNA = []
    originalTRN = []
    #lists for the powersets of the TRNA and then a list containing just everything set but in string form.
    def __init__(self, header, sequence):
        '''
        In this __init__ method the sequence is stripped from any '.' and '-' that be within the FastA files
        It also creates the powersets that are to be stored in by lists within uniqueTRNA. This way the rest of the class can call upon the lists without emptying it
        '''
        self.sequence = sequence
        self.sequence = self.sequence.replace('.','')
        self.sequence = self.sequence.replace('-','')
        self.header = header
        uniqueTRNA.powerset = set()
        uniqueTRNA.originalTRN.append(sequence)
        for i in range(len(sequence)):
            for j in range(i, len(sequence)):
                setseq=(sequence[i:j+1])
                self.powerset.add(setseq)
                self.originalTRN.append(setseq)
        self.TRNA.append(self.powerset)
    def setgenerator(self):
        '''
        Will create and store different sets for the uniques and essential sequences.
        The uniques are found by repeating through the powersets and only adding the unqiue sets that appear once, to the unique set.
        The essentials are found by repeatimg each item within the unique set and seaching if appear in all set or set's substring
        '''
        nonessentials = set()
        amountpowerset = list(range(0,len(uniqueTRNA.TRNA)))
        uniqset = uniqueTRNA.TRNA[0]
        for n in range(1,len(amountpowerset)):
            uniqset = uniqset^uniqueTRNA.TRNA[n]

        preessentials = []
        for cont in uniqset:
            if uniqueTRNA.originalTRN.count(cont) == 1:
                preessentials.append(cont)

        subessentials = []

        for item in preessentials:
            for item2 in preessentials:
                if item in item2 and item != item2:
                    subessentials.append(item2)

        essentials = set(preessentials) ^ set(subessentials)
        finalessentials=list(essentials)

        return finalessentials

def main(inCL=None):
    '''
    This main function prints each header and sequence then the uniqueessentials sequences to organizes the uniqueessentials data.
    The sequences are organized by their header alphabetically.
    Example.
    ABCD
    ..CD
    '.' are placed where no uniqueessential sequence is present. When a uniqueessential is present then it's printed.
    '''
    myReader = FastAreader()
    costcolist = []
    for header, sequence in myReader.readFasta():
        sample = uniqueTRNA(header,sequence)

        costcolist.append(sample)
    sorted(costcolist,key = lambda x: x.header)
    uniqueessentials = sample.setgenerator()

    for sample in costcolist:
        print(sample.header)
        print(sample.sequence)
        for step in uniqueessentials:
            dots = sample.sequence.find(step)
            if dots > -1:
                print(dots*'.'+step)



if __name__ == "__main__":

    main()


# In[ ]:




