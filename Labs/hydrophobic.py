# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A

from sequenceAnalysis import FastAreader
from sequenceAnalysis import CommandLine
from sequenceAnalysis import OrfFinder
from sequenceAnalysis import Hydrophobicity

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure


def main(inCL=None):

    if inCL is None:
        myCommandLine = CommandLine() #['tass2.fa','tass2ORFdata-ATG-100.txt','--longestGene']
    else :
        myCommandLine = CommandLine(inCL)

    longest = myCommandLine.args.longestGene # is True if only the longest Gene is desired
    start = myCommandLine.args.start # is a list of start codons
    minGene = myCommandLine.args.minGene # is the minimum Gene length to include
    stop = myCommandLine.args.stop # optional stop arguments

    # setup our orfReader
    orfReader = FastAreader(myCommandLine.args.inFile)  # runs FastAreader on command line infile
    for head, seq in orfReader.readFasta():
        orf = OrfFinder(seq, head, start, stop, minGene, longest)
        orf.orf()


    hydroReader = FastAreader(myCommandLine.args.inFile)
    with open(myCommandLine.args.outFile, 'w') as outFile:
        for head, seq in hydroReader.readFasta():
            for item in orf.orfData:
                if item[0] == 1:
                    rnaSeq = seq.replace('T', 'U')
                    beginning = int(item[1])
                    end = int(item[2])
                    hydro = Hydrophobicity(rnaSeq[beginning - 1:end], head, item[0], item[1], item[2], item[3])
                    hydro.out(outFile)
                if item[0] == 2:
                    rnaSeq = seq.replace('T', 'U')
                    beginning = int(item[1])
                    end = int(item[2])
                    hydro = Hydrophobicity(rnaSeq[beginning - 1:end], head, item[0], item[1], item[2], item[3])
                    hydro.out(outFile)
                if item[0] == 3:
                    rnaSeq = seq.replace('T', 'U')
                    beginning = int(item[1])
                    end = int(item[2])
                    hydro = Hydrophobicity(rnaSeq[beginning - 1:end], head, item[0], item[1], item[2], item[3])
                    hydro.out(outFile)

    #https://gist.github.com/marcelm/c0cbb0b6ee44b471b910
    with PdfPages('multi.pdf') as pages:
        for head, seq in hydroReader.readFasta():
            for item in orf.orfData:
                if item[0] == 1:
                    rnaSeq = seq.replace('T', 'U')
                    beginning = int(item[1])
                    end = int(item[2])
                    hydro = Hydrophobicity(rnaSeq[beginning - 1:end], head, item[0], item[1], item[2], item[3])

                    L = []
                    for item in hydro.aaToValue:
                        if item == '-STOP-':
                            hydro.aaToValue.remove('-STOP-')
                        else:
                            continue

                    if len(hydro.aaToValue) > 1000: # if range greater than 1000 choose a window range (20?)
                        x = np.convolve(hydro.aaToValue, np.ones(45), 'valid')/45

                    elif len(hydro.aaToValue) > 500: # if range greater than 500 (but less than 100)
                        x = np.convolve(hydro.aaToValue, np.ones(30), 'valid')/30

                    elif len(hydro.aaToValue) > 100: # if range greater than 100 (but less than 500)
                        x = np.convolve(hydro.aaToValue, np.ones(15), 'valid')/15

                    else:                           # if range less than 100
                        x = np.convolve(hydro.aaToValue, np.ones(5), 'valid')/5

                    #y = np.array(L)
                    fig = Figure()
                    ax = fig.gca()
                    ax.plot(x)
                    ax.set_ylabel('Hydrophobicity Average')
                    ax.set_xlabel('Amino Acid Position')
                    ax.set_title('Frame 1')

                    canvas = FigureCanvasPdf(fig)
                    canvas.print_figure(pages)

                    L.clear()

                if item[0] == 2:
                    rnaSeq = seq.replace('T', 'U')
                    beginning = int(item[1])
                    end = int(item[2])
                    hydro = Hydrophobicity(rnaSeq[beginning - 1:end], head, item[0], item[1], item[2], item[3])

                    for item in hydro.aaToValue:
                        if item == '-STOP-':
                            hydro.aaToValue.remove('-STOP-')
                        else:
                            continue

                    if len(hydro.aaToValue) > 1000: # if range greater than 1000 choose a window range (20?)
                        x = np.convolve(hydro.aaToValue, np.ones(45), 'valid')/45

                    elif len(hydro.aaToValue) > 500: # if range greater than 500 (but less than 100)
                        x = np.convolve(hydro.aaToValue, np.ones(30), 'valid')/30

                    elif len(hydro.aaToValue) > 100: # if range greater than 100 (but less than 500)
                        x = np.convolve(hydro.aaToValue, np.ones(15), 'valid')/15

                    else:                           # if range less than 100
                        x = np.convolve(hydro.aaToValue, np.ones(5), 'valid')/5

                    #y = np.array(L)
                    fig = Figure()
                    ax = fig.gca()
                    ax.plot(x)
                    ax.set_ylabel('Hydrophobicity Average')
                    ax.set_xlabel('Amino Acid Position')
                    ax.set_title('Frame 2')

                    canvas = FigureCanvasPdf(fig)
                    canvas.print_figure(pages)

                    L.clear()

                if item[0] == 3:
                    rnaSeq = seq.replace('T', 'U')
                    beginning = int(item[1])
                    end = int(item[2])
                    hydro = Hydrophobicity(rnaSeq[beginning - 1:end], head, item[0], item[1], item[2], item[3])

                    for item in hydro.aaToValue:
                        if item == '-STOP-':
                            hydro.aaToValue.remove('-STOP-')
                        else:
                            continue

                    if len(hydro.aaToValue) > 1000: # if range greater than 1000 choose a window range (20?)
                        x = np.convolve(hydro.aaToValue, np.ones(45), 'valid')/45


                    elif len(hydro.aaToValue) > 500: # if range greater than 500 (but less than 100)
                        x = np.convolve(hydro.aaToValue, np.ones(30), 'valid')/30

                    elif len(hydro.aaToValue) > 100: # if range greater than 100 (but less than 500)
                        x = np.convolve(hydro.aaToValue, np.ones(15), 'valid')/15

                    else:                           # if range less than 100
                        x = np.convolve(hydro.aaToValue, np.ones(5), 'valid')/5

                    #y = np.array(x)
                    fig = Figure()
                    ax = fig.gca()
                    ax.plot(x)
                    ax.set_ylabel('Hydrophobicity Average')
                    ax.set_xlabel('Amino Acid Position')
                    ax.set_title('Frame 3')

                    canvas = FigureCanvasPdf(fig)
                    canvas.print_figure(pages)

                    L.clear()


if __name__ == "__main__":
    main()
