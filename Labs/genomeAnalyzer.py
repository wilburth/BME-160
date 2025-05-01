# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A


from sequenceAnalysis import NucParams          # selective import !
from sequenceAnalysis import FastAreader        # selective import !

def main ():

    myReader = FastAreader('testGenome.fa')                    # make sure to change this to use stdin
    myNuc = NucParams()                         # here is where we rename our selective import
    for head, seq in myReader.readFasta():
        myNuc.addSequence(seq)                  # call addSequence method from class NucParams on 'seq'

    # find sequence length
    # use aaComp to count amino acids, then multiply by 3
    # return count divided by 1,000,000 -->  per million bases
    aaComp = myNuc.aaComposition()   # set aaComp to be our amino acid dictionary from class NucParams
    nucCount = 0                     # initialize nucleotide count at 0
    for key in aaComp:               # iterate across all keys in the amino acid dictionary
        nucCount += (aaComp[key]*3)  # for each key, multiply its matching value times 3, summate all values, set nucCount equal to this summated value (for each amino acid to occur / be coded, 3 nucleotides must exist !)

    print('Sequence length = {:0.2f} Mb'.format(nucCount/1000000), '\n')

    # find GC content
    # in nucleotide dictionary --> we have {key:value} --> sum value for G + C , then divide by total sum of values
    nucComp = myNuc.nucComposition() # set nucComp to nucleotide dictionary from class NucParams

    gcCount = nucComp['G'] + nucComp['C']
    totalCount = sum(nucComp.values())

    print('GC Content = {0:0.1f}%'.format((gcCount/totalCount) * 100), '\n')

    # sort codons in alpha order, by Amino Acid
    # https://docs.python.org/3/howto/sorting.html
    # https://www.w3schools.com/python/python_lambda.asp
    # https://stackoverflow.com/questions/17341421/list-with-multiple-values-per-key
    # nucs = sorted(myNuc.rnaCodonTable.items(), Lambda functions can take any number of arguments sort: by value, then by key)
    # sort by value, then by letter --> will print '-' first, then sort in alphabetical (for amino acids)
    # have to sort rnaCodonTable ITEMS (contains the key-value pairs of the dictionary, as tuples in a list)
    nucs = sorted(myNuc.rnaCodonTable.items(), key = lambda values: (values[1], values[0]))

    #for item in nucs:
    #    print(item, '\n')

    # calculate relative codon usage for each codon and print
    codonComp = myNuc.codonComposition()        # set codonComp to be our codon dictionary from class NucParams
    for nucI in nucs:       # for all item in nucs (sorted list now working !)

        nuc = nucI[0] # first position for each item (example: UAA)
        aa = nucI[1]  # second position for each item (example: - )

        # amount of times the codon occurs divided by amount of times amino acid occurs (since some amino acids are coded by multiple codons) (example: 32.6) in percent format
        val = codonComp[nuc] / aaComp[aa]

        print ('{:s} : {:s} {:5.1f} ({:6d})'.format(nuc, aa, val*100, codonComp[nuc]))

if __name__ == "__main__":
    main()
