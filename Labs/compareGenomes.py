# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A



from sequenceAnalysis import NucParams
from sequenceAnalysis import FastAreader


def main():
    """

    between TWO (2) genomes:
        compare GC content
        compare aaComposition
        compare relative codon bias

    halophile genome + hyperthermophile genome to compare

    """

    myReader1 = FastAreader('testGenome.fa')
    myNuc1 = NucParams()
    for head1, seq1 in myReader1.readFasta():
        myNuc1.addSequence(seq1)

    myReader2 = FastAreader('haloVolc1_1-genes.fa')
    myNuc2 = NucParams()
    for head2, seq2 in myReader2.readFasta():
        myNuc2.addSequence(seq2)

    """"COMPARE NUCLEOTIDE COUNTS"""
    aaComp1 = myNuc1.aaComposition()
    nucCount1 = 0
    for key in aaComp1:
        nucCount1 += (aaComp1[key]*3)

    aaComp2 = myNuc2.aaComposition()
    nucCount2 = 0
    for key in aaComp2:
        nucCount2 += (aaComp2[key]*3)

    print('')
    print('Sequence length of (first genome) | (second genome):      {:0.2f} Mb | {:0.2f} Mb'.format((nucCount1/1000000), (nucCount2/1000000)))

    """COMPARE GC CONTENT HERE"""
    nucComp1 = myNuc1.nucComposition()
    nucComp2 = myNuc2.nucComposition()

    gcCount1 = nucComp1['G'] + nucComp1['C']
    totalCount1 = sum(nucComp1.values())
    gcCount2 = nucComp2['G'] + nucComp2['C']
    totalCount2 = sum(nucComp2.values())

    gcContent1 = (gcCount1/totalCount1) * 100
    gcContent2 = (gcCount2/totalCount2) * 100

    print('')
    print('GC Content of (first genome) | (second genome):      {:0.1f}% | {:0.1f}%'.format((gcContent1), gcContent2))


    """COMAPRE RELATIVE CODON BIAS"""

    print('')
    print('Codon usage in % format \nCodon : aa : (first genome) <=> (second genome)')
    nucs = sorted(myNuc1.rnaCodonTable.items(), key = lambda values: (values[1], values[0]))
    codonComp1 = myNuc1.codonComposition()
    codonComp2 = myNuc2.codonComposition()
    for nucI in nucs:

        nuc = nucI[0]
        aa = nucI[1]

        val1 = codonComp1[nuc] / aaComp1[aa]
        val2 = codonComp2[nuc] / aaComp2[aa]

        if val1 > val2:
            print('{:s} : {:s} : {:5.1f}  >  {:5.1f}'.format(nuc, aa, val1*100, val2*100))
        elif val1 < val2:
            print('{:s} : {:s} : {:5.1f}  <  {:5.1f}'.format(nuc, aa, val1*100, val2*100))
        else:
            print('{:s} : {:s} : {:5.1f}  =  {:5.1f}'.format(nuc, aa, val1*100, val2*100))

    print('')
    print('Amino Acid count \nCodon : aa : (first genome) <=> (second genome)')
    for nucI in nucs:

        nuc = nucI[0]
        aa = nucI[1]

        if codonComp1[nuc] > codonComp2[nuc]:
            print('{:s} : {:s} : ({:6d})  >  ({:6d})'.format(nuc, aa, codonComp1[nuc], codonComp2[nuc]))
        elif codonComp1[nuc] < codonComp2[nuc]:
            print('{:s} : {:s} : ({:6d})  <  ({:6d})'.format(nuc, aa, codonComp1[nuc], codonComp2[nuc]))
        else:
            print('{:s} : {:s} : ({:6d})  =  ({:6d})'.format(nuc, aa, codonComp1[nuc], codonComp2[nuc]))





if __name__ == "__main__":
    main()