# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A

from sequenceAnalysis import FastAreader


# will need methods to create power set, find the uniques
# also want method to clean the sequence
# also want method for printing output --> only use main to call methods from class FindUnique
#

##### following rough design plan provided in lab06 by Dr. B ######
#
# compute the set of all substrings from each tRNA sequence.
# [powerset] I will refer to a set of substrings as a tRNAset.
#
# for each tRNAset, compute the union of all other tRNA sets, and then remove that union from the current tRNAset.
# Notice that this union operation finds all of the substrings from all other tRNA.
# If any of those are present in your current tRNA, then they are not unique ! [ Unique]
#
# for each unique tRNAset, it now contains the truly unique ones, along with any extensions of that subsequence.
# If, for example, it was found that G only occurred in a single tRNA, then adding an A onto that G must also be unique
# because it has a G in it. We only want the minimal form.. G. [ Essential]
#
# Remove spaces from the header line before sorting and printing. This will make your output a little prettier.
#
# make sure to remove any alignment characters from your initial sequences.
# These characters are periods(.), underscores(\_), or dashes(-) .
# You will find many new characters in the sequence other than {ACGU} - leave these in place, they are modified bases.


class FindUnique:
    """
    FindUnique docstring: Class FindUnique, which is within program findUnique.py, is to be used to find
    unique and essential tRna while reading in from a .fa (fasta) file.

    """
    def __init__(self):
        """
        __init__() docstring: __init__() declares / instantiates all necessary objects to be used by the other methods
        within class FindUnique. Additionally, __init__() handles the reading of the fasta file from stdin.

        """
        # need: powerset list, uniques list --> try out dictionary filled with header + sequence from FastAFIle
        self.powSetList = []
        self.uniquesList = []
        self.headSeqDict = {}

        myReader = FastAreader('bos-tRNA-7.fa') # instantiate fastareader # 'bos-tRNA-7.fa' remove when you want to run with stdin

        pos = 0 # starting item position (within dictionary)
        for head, seq in myReader.readFasta():
            cleanSeq = self.clean(seq) # call method cleanSeq to clean the sequence from fastaFile
            self.headSeqDict[pos] = [head, cleanSeq]    # position in headSeq dictionary (starting at pos=0) --> at each position, fill it with list containing the header + sequence from that fasta Line
            self.powSetList.append(self.powerSet(cleanSeq))     # create power sets by calling method 'powerSet()' on 'cleanSeq' --> append it to powSetList
            pos += 1 # count up in position one

        #print(self.headSeqDict) # working
        #print(self.powSetList) # working

    def powerSet(self, seq):
        """
        powerSet() docstring: powerSet() is called in the __init__() method on the cleaned sequence from the Fasta File,
        and is used to create power sets of all the characters in the sequence.

        """
        pow = set() # need to make use of sets

        for i in range(len(seq)): # for i in range of length of seq
            length_var2 = len(seq)  # length_var2 (also set equal to length of seq)

            while length_var2 > i:  # while length_var2 is greater than position 'i'
                pow.add(seq[i: length_var2])    # add to our 'pow' set the sequence, cut from position 'i' up to length_var2
                length_var2 -= 1        # count down in length by 1

        return pow  # return 'pow' of type set

    def unique(self):
        """
        unique() docstring: unique() makes multiple copies of the power sets in order to distinguish duplicates and
        remove them to ultimately obtain all of the uniques that occur in the tRna sequence

        """
        for sets in self.powSetList:
            union = set() # object 'union' of type set
            listCopy = self.powSetList.copy() # copyList set to a copy of powSetList (which has been appended with all power sets)
            setCopy = sets.copy()   # copySet set to a copy of a set in powSetList
            listCopy.remove(setCopy) # remove from copyList whatever is also in copySet

            #print(copySet)

            for powSet in listCopy: # for each power set in copyList
                union = union.union(powSet) # use object 'union', call function union() with parameter powSet --> https://www.w3schools.com/python/ref_set_union.asp
            setCopy.difference_update(union) # compute difference on copySet with parameter as object 'union' --> https://www.w3schools.com/python/ref_set_difference_update.asp

            copySet_obj2 = setCopy.copy() # object 'copySet_obj2' which is a copy of a copy of the sets in powSetList with the difference_update method applied to it
            for tRna in setCopy:
                uniques = setCopy.copy()
                uniques.remove(tRna) # pull out tRna from uniques

                for tRna_2 in uniques: # if tRna #2 is in uniques
                    if tRna in tRna_2 and len(tRna) < len(tRna_2): # and if tRna #1 is in tRna #2, along with tRna #2 being longer
                        copySet_obj2.discard(tRna_2) # remove it from our copySet 2, we keep the essentials--> https://www.w3schools.com/python/ref_set_discard.asp

            self.uniquesList.append(copySet_obj2) # append copySet_obj2 to uniquesList

    def clean(self, seq):
        """
        clean() docstring: a quick solution to cleaning the string of specific unwanted characters, can be used anywhere

        """
        clean = seq.replace('_', '').replace('-','').replace('.', '') # "They are not part of the sequence and must be removed. [-\_\.] are alignment characters."
        return clean

    def out(self):
        """
        out() docstring: out() parses through the headSeqDict, which was created in the __init()__ method, and uses it
        to print the header and sequence for each tRna entry in an input fasta file, followed by a formatted printing
        of the uniques
        """



        #yuh = sorted(self.headSeqDict.values())
        #print(yuh)

        #for index, object in enumerate(yuh):
            #print(yuh[index][0])
            #print(yuh[index][1])

            #for pos in range(len(yuh[index][1])):
                #for string in self.uniquesList[index]:
                    #if string == yuh[index][1][pos: pos + len(string)]:
                        #print(('.' * pos) + string)

        ### following block of code is updated to print fasta input file in alphabetical order ###
        sortedDict = dict(sorted(self.headSeqDict.items(), key=lambda i:i[1]))
        for index in sortedDict:
            print(sortedDict[index][0], '\n', sortedDict[index][1])

            for pos in range(len(sortedDict[index][1])):
                for string in self.uniquesList[index]:
                    if string == sortedDict[index][1][pos: pos + len(string)]:
                        print(('.' * pos) + string)

        #for item in self.headSeqDict: # for each item in our dictionary
            #print(item)
            #print(self.headSeqDict[item][0], '\n', self.headSeqDict[item][1]) # prints the header, followed by a new line with the sequence on it for each item

            #for pos in range(len(self.headSeqDict[item][1])): # for each position in the sequence (of each item)
                #for string in self.uniquesList[item]: # for each string in our uniquesList (at position 'item')
                    #if string == self.headSeqDict[item][1][pos: pos + len(string)]: # if the string matches the sequence which is cut from position x up to position y (which will be the length of the string)
                        #print(('.' * pos) + string) # print the correct amount of periods followed by the string


def main():
    '''
    main() docstring: main() instantiates object 'tRNA' to be of class FindUnique, then calls methods
    unique() and out() on object 'tRNA'.

    '''

    tRNA = FindUnique()
    tRNA.unique()
    tRNA.out()


if __name__ == "__main__":
    main()